#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include "rtklib.h"
#include "vbs_core.h"
#include "vbs_merge.h"
#include "tcp_io.h"
#include "net_interp.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

typedef struct {
    const char *src_host; int src_port; int dual; const char *b_host; int b_port; int epoch_window_ms;
    const char *eph_host; int eph_port; const char *ssr_host; int ssr_port; int out_port;
    int use_llh_offset; double dn, de, du, vbs_lat, vbs_lon, vbs_alt;
    int apply_trop, iono_interp, msm_out_type, eph_out_enable, ssr_passthrough;
    int trace_level; const char *trace_file;
} config_t;

typedef struct { tcp_client_t tcp; rtcm_t rtcm; const char *host; int port; const char *tag; long last_conn_try_ms; } channel_t;
static volatile int g_stop = 0; static void on_sigint(int sig){(void)sig; g_stop=1;}
static int is_ssr_msg(int t){return((t>=1057&&t<=1068)||(t>=1240&&t<=1254)||(t>=1258&&t<=1270)||t==11||t==12||t==13||t==14);}
static int frame_mtype(const uint8_t *b,int l){return(l<6||b[0]!=0xD3)?-1:(int)getbitu(b,24,12);}
static int msm_type_for_sys(int s,int st){switch(s){case SYS_GPS:return 1070+st;case SYS_GLO:return 1080+st;case SYS_GAL:return 1090+st;case SYS_SBS:return 1100+st;case SYS_QZS:return 1110+st;case SYS_CMP:return 1120+st;case SYS_IRN:return 1130+st;default:return 0;}}
static int eph_msg_for_sys(int s){switch(s){case SYS_GPS:return 1019;case SYS_GLO:return 1020;case SYS_QZS:return 1044;case SYS_GAL:return 1046;case SYS_CMP:return 1042;case SYS_IRN:return 1041;default:return 0;}}
static long now_ms_mono(void);
static void channel_init(channel_t*c,const char*h,int p,const char*t){c->tcp.fd=-1;c->host=h;c->port=p;c->tag=t;c->last_conn_try_ms=0;init_rtcm(&c->rtcm);c->rtcm.time.time=time(NULL);c->rtcm.time.sec=0.0;}
static void channel_close(channel_t*c){tcpc_close(&c->tcp);free_rtcm(&c->rtcm);}
static void channel_try_connect(channel_t*c){long now=now_ms_mono();if(c->tcp.fd>=0)return;if(c->last_conn_try_ms&&now-c->last_conn_try_ms<1000)return;c->last_conn_try_ms=now;if(tcpc_connect(&c->tcp,c->host,c->port)==0)fprintf(stderr,"[%s] connected %s:%d\n",c->tag,c->host,c->port);}
static long now_ms_mono(void)
{
#ifdef _WIN32
    return (long)GetTickCount64();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long)(ts.tv_sec * 1000 + ts.tv_nsec / 1000000);
#endif
}

#ifdef _WIN32
static void sleep_ms(int ms) { Sleep((DWORD)ms); }
#else
static void sleep_ms(int ms)
{
    struct timespec ts = {0, (long)ms * 1000000L};
    nanosleep(&ts, NULL);
}
#endif

static void setup_signals(void)
{
    signal(SIGINT, on_sigint);
    signal(SIGTERM, on_sigint);
#ifndef _WIN32
    signal(SIGPIPE, SIG_IGN);
#endif
}
static void mirror_eph(rtcm_t*out,const rtcm_t*in){int sat=in->ephsat,set=in->ephset,prn=0,sys=satsys(sat,&prn);if(sys==SYS_GLO)out->nav.geph[prn-1]=in->nav.geph[prn-1];else if(sys!=0)out->nav.eph[sat-1+MAXSAT*set]=in->nav.eph[sat-1+MAXSAT*set];out->ephsat=sat;out->ephset=set;}
static void mirror_glo_fcn(rtcm_t*out,const rtcm_t*in){memcpy(out->nav.glo_fcn,in->nav.glo_fcn,sizeof(in->nav.glo_fcn));}
static void copy_obs(obsd_t*d,int*n,const obsd_t*s,int m){int k=m>MAXOBS?MAXOBS:m;memcpy(d,s,k*sizeof(obsd_t));*n=k;}
static void emit_station(rtcm_t*out,tcp_server_t*srv){if(gen_rtcm3(out,1006,0,0))tcps_broadcast(srv,out->buff,out->nbyte);}
static void emit_eph(rtcm_t*out,tcp_server_t*srv,int sat){int sys=satsys(sat,NULL),msg=eph_msg_for_sys(sys);if(!msg)return; if(sys==SYS_GAL){if(gen_rtcm3(out,1046,0,0))tcps_broadcast(srv,out->buff,out->nbyte);if(gen_rtcm3(out,1045,0,0))tcps_broadcast(srv,out->buff,out->nbyte);}else if(gen_rtcm3(out,msg,0,0))tcps_broadcast(srv,out->buff,out->nbyte);}
static void emit_msm(rtcm_t*out,tcp_server_t*srv,int st){obsd_t*data=out->obs.data,master[MAXOBS],sys_buf[MAXOBS];int nobs=out->obs.n,mn=nobs<MAXOBS?nobs:MAXOBS;static const int syses[]={SYS_GPS,SYS_GLO,SYS_GAL,SYS_CMP,SYS_QZS,SYS_IRN,SYS_SBS};if(nobs<=0)return;memcpy(master,data,mn*sizeof(obsd_t));int total=0;for(size_t k=0;k<sizeof(syses)/sizeof(syses[0]);k++){int sys=syses[k],nsat=0,nsig=0,mask[MAXCODE]={0};if(!msm_type_for_sys(sys,st))continue;for(int i=0;i<mn;i++)if(satsys(master[i].sat,NULL)==sys){nsat++;for(int j=0;j<NFREQ+NEXOBS;j++){int code=master[i].code[j];if(code&&!mask[code-1]){mask[code-1]=1;nsig++;}}}if(nsat>0&&nsig>0&&nsig<=64)total+=(nsat-1)/(64/nsig)+1;}if(total<=0)return;int sent=0;for(size_t k=0;k<sizeof(syses)/sizeof(syses[0]);k++){int sys=syses[k],msg=msm_type_for_sys(sys,st),nsat=0,nsig=0,mask[MAXCODE]={0};if(!msg)continue;for(int i=0;i<mn;i++)if(satsys(master[i].sat,NULL)==sys){sys_buf[nsat++]=master[i];for(int j=0;j<NFREQ+NEXOBS;j++){int code=master[i].code[j];if(code&&!mask[code-1]){mask[code-1]=1;nsig++;}}}if(nsat<=0||nsig<=0||nsig>64)continue;int ns=64/nsig,nmsg=(nsat-1)/ns+1,off=0;out->obs.data=sys_buf;for(int m=0;m<nmsg;m++){int chunk=(nsat-off)<ns?(nsat-off):ns;if(off>0)memmove(sys_buf,sys_buf+off,chunk*sizeof(obsd_t));out->obs.n=chunk;if(gen_rtcm3(out,msg,0,(sent<total-1)?1:0))tcps_broadcast(srv,out->buff,out->nbyte);off+=chunk;sent++;}out->obs.data=data;out->obs.n=nobs;}}
static void usage(const char*p){fprintf(stderr,"Usage: %s [options]\n--dual --b-port PORT --iono --trop --offset N E U | --vbs LAT LON ALT\n",p);}
static void run(const config_t*cfg){if(cfg->trace_file)traceopen(cfg->trace_file);tracelevel(cfg->trace_level);channel_t chA,chB,chEph,chSsr;channel_init(&chA,cfg->src_host,cfg->src_port,"A");if(cfg->dual)channel_init(&chB,cfg->b_host,cfg->b_port,"B");if(cfg->eph_port)channel_init(&chEph,cfg->eph_host,cfg->eph_port,"eph");if(cfg->ssr_port)channel_init(&chSsr,cfg->ssr_host,cfg->ssr_port,"ssr");rtcm_t rout;init_rtcm(&rout);rout.time.time=time(NULL);vbs_ctx_t vbs;double vbs_llh[3]={cfg->vbs_lat,cfg->vbs_lon,cfg->vbs_alt};vbs_init(&vbs,cfg->use_llh_offset,cfg->dn,cfg->de,cfg->du,cfg->use_llh_offset?NULL:vbs_llh,cfg->apply_trop);merger_t mrg;merger_init(&mrg,cfg->epoch_window_ms);net_interp_t neti;net_interp_init(&neti,cfg->iono_interp);tcp_server_t*srv=tcps_create(cfg->out_port,16);if(!srv){fprintf(stderr,"[error] cannot bind output port %d\n",cfg->out_port);return;}int subtype=cfg->msm_out_type?cfg->msm_out_type:7;long last_stat_t=(long)time(NULL),ssr_fwd=0;obsd_t last_a[MAXOBS],last_b[MAXOBS];int last_a_n=0,last_b_n=0;gtime_t last_a_t={0},last_b_t={0};int side_staid[2]={0,0};uint8_t buf[4096];while(!g_stop){channel_try_connect(&chA);if(cfg->dual)channel_try_connect(&chB);if(cfg->eph_port)channel_try_connect(&chEph);if(cfg->ssr_port)channel_try_connect(&chSsr);int any_data=0;channel_t*chs[4]={&chA,NULL,NULL,NULL};int nch=1;if(cfg->dual)chs[nch++]=&chB;if(cfg->eph_port)chs[nch++]=&chEph;if(cfg->ssr_port)chs[nch++]=&chSsr;for(int ci=0;ci<nch;ci++){channel_t*ch=chs[ci];if(ch->tcp.fd<0)continue;for(;;){int r=tcpc_read(&ch->tcp,buf,sizeof(buf));if(r==TCPC_AGAIN)break;if(r==0||r<0){if(r!=0||ch->tcp.fd>=0)fprintf(stderr,"[%s] disconnected\n",ch->tag);tcpc_close(&ch->tcp);break;}any_data=1;for(int i=0;i<r;i++){int ret=input_rtcm3(&ch->rtcm,buf[i]);if(!ret)continue;int mtype=frame_mtype(ch->rtcm.buff,ch->rtcm.len+3);if(is_ssr_msg(mtype)){if(cfg->ssr_passthrough){tcps_broadcast(srv,ch->rtcm.buff,ch->rtcm.len+3);ssr_fwd++;}continue;}rout.time=ch->rtcm.time;rout.staid=ch->rtcm.staid;if(ret==1){mirror_glo_fcn(&rout,&ch->rtcm);if(cfg->dual&&(ch==&chA||ch==&chB)){int side=(ch==&chA)?0:1;side_staid[side]=ch->rtcm.staid;merger_stash(&mrg,side,ch->rtcm.obs.data,ch->rtcm.obs.n,ch->rtcm.time,now_ms_mono());if(!side){copy_obs(last_a,&last_a_n,ch->rtcm.obs.data,ch->rtcm.obs.n);last_a_t=ch->rtcm.time;}else{copy_obs(last_b,&last_b_n,ch->rtcm.obs.data,ch->rtcm.obs.n);last_b_t=ch->rtcm.time;}}else if(ch==&chA){int n=ch->rtcm.obs.n>MAXOBS?MAXOBS:ch->rtcm.obs.n;memcpy(rout.obs.data,ch->rtcm.obs.data,n*sizeof(obsd_t));rout.obs.n=n;vbs_correct_obs(&vbs,&rout,NULL);if(rout.obs.n>0)emit_msm(&rout,srv,subtype);}}else if(ret==2){mirror_eph(&rout,&ch->rtcm);if(cfg->eph_out_enable)emit_eph(&rout,srv,ch->rtcm.ephsat);}else if(ret==5){int side=(ch==&chA)?0:((cfg->dual&&ch==&chB)?1:-1);if(side>=0){side_staid[side]=ch->rtcm.staid;vbs_update_base(&vbs,side,&ch->rtcm.sta);if(norm(vbs.vbs_ecef,3)>=1e3){rout.sta=ch->rtcm.sta;vbs_rewrite_station(&vbs,&rout);emit_station(&rout,srv);}}}}}}if(cfg->dual){obsd_t raw_a[MAXOBS],raw_b[MAXOBS];int na=0,nb=0,have_a=0,have_b=0;gtime_t mt;double sat_extra_m[MAXSAT]={0};if(merger_poll_pair(&mrg,now_ms_mono(),raw_a,&na,&have_a,raw_b,&nb,&have_b,&mt)&&(have_a||have_b)){if(cfg->iono_interp&&have_a&&have_b)net_interp_build_dual(&neti,&vbs,&rout,raw_a,na,raw_b,nb,sat_extra_m);static unsigned char side_a[MAXOBS],side_b[MAXOBS];if(have_a){memset(side_a,0,sizeof(side_a));memcpy(rout.obs.data,raw_a,na*sizeof(obsd_t));rout.obs.n=na;rout.time=mt;vbs_correct_obs_ex(&vbs,&rout,side_a,cfg->iono_interp?sat_extra_m:NULL);na=rout.obs.n;memcpy(raw_a,rout.obs.data,na*sizeof(obsd_t));if(na==0)have_a=0;}if(have_b){memset(side_b,1,sizeof(side_b));memcpy(rout.obs.data,raw_b,nb*sizeof(obsd_t));rout.obs.n=nb;rout.time=mt;vbs_correct_obs_ex(&vbs,&rout,side_b,cfg->iono_interp?sat_extra_m:NULL);nb=rout.obs.n;memcpy(raw_b,rout.obs.data,nb*sizeof(obsd_t));if(nb==0)have_b=0;}if(have_a||have_b){if(have_a&&side_staid[0])rout.staid=side_staid[0];else if(have_b&&side_staid[1])rout.staid=side_staid[1];int nm=merger_align_and_merge(&mrg,raw_a,na,have_a,raw_b,nb,have_b,rout.obs.data);rout.obs.n=nm;rout.time=mt;if(nm>0)emit_msm(&rout,srv,subtype);}}}if(!any_data){sleep_ms(5);}long now_s=(long)time(NULL);if(now_s-last_stat_t>=15){if(cfg->dual)fprintf(stderr,"[stat] 1005/6=%ld MSM=%ld | obs in=%ld out=%ld | sat corr=%ld no_eph=%ld net=%ld | epochs merged=%ld A-only=%ld B-only=%ld | CNR both_sat=%ld chose_A=%ld chose_B=%ld | align init=%ld reset(slip=%ld code=%ld gap=%ld resid=%ld) | Bnorm=%ld Bunal=%ld switch_blk=%ld | iono epochs=%ld pairs=%ld interp=%ld rej=%ld | SSR fwd=%ld | clients=%d\n",vbs.n_frames_1005_6,vbs.n_frames_msm,vbs.n_obs_in,vbs.n_obs_out,vbs.n_sat_corrected,vbs.n_sat_no_eph,vbs.n_sat_netcorr,mrg.n_merged,mrg.n_a_only,mrg.n_b_only,mrg.n_both,mrg.n_chose_a,mrg.n_chose_b,mrg.n_align_init_ok,mrg.n_align_reset_slip,mrg.n_align_reset_code,mrg.n_align_reset_gap,mrg.n_align_reset_resid,mrg.n_b_norm_applied,mrg.n_b_unaligned_lli,mrg.n_switch_blocked,neti.n_epochs_used,neti.n_sat_pairs,neti.n_sat_interp,neti.n_sat_reject,ssr_fwd,tcps_num_clients(srv));else fprintf(stderr,"[stat] 1005/6=%ld MSM=%ld | obs in=%ld out=%ld | sat corr=%ld no_eph=%ld | SSR fwd=%ld | clients=%d\n",vbs.n_frames_1005_6,vbs.n_frames_msm,vbs.n_obs_in,vbs.n_obs_out,vbs.n_sat_corrected,vbs.n_sat_no_eph,ssr_fwd,tcps_num_clients(srv));last_stat_t=now_s;}}channel_close(&chA);if(cfg->dual)channel_close(&chB);if(cfg->eph_port)channel_close(&chEph);if(cfg->ssr_port)channel_close(&chSsr);tcps_destroy(srv);free_rtcm(&rout);}
int main(int argc,char**argv){config_t cfg={.src_host="127.0.0.1",.src_port=50001,.dual=0,.b_host="127.0.0.1",.b_port=50003,.epoch_window_ms=40,.eph_host="127.0.0.1",.eph_port=0,.ssr_host="127.0.0.1",.ssr_port=0,.out_port=50002,.use_llh_offset=1,.dn=0,.de=0,.du=0,.vbs_lat=0,.vbs_lon=0,.vbs_alt=0,.apply_trop=0,.iono_interp=0,.msm_out_type=7,.eph_out_enable=1,.ssr_passthrough=1,.trace_level=0,.trace_file=NULL};for(int i=1;i<argc;i++){if(!strcmp(argv[i],"--src-host")&&i+1<argc)cfg.src_host=argv[++i];else if(!strcmp(argv[i],"--src-port")&&i+1<argc)cfg.src_port=atoi(argv[++i]);else if(!strcmp(argv[i],"--dual"))cfg.dual=1;else if(!strcmp(argv[i],"--b-host")&&i+1<argc)cfg.b_host=argv[++i];else if(!strcmp(argv[i],"--b-port")&&i+1<argc)cfg.b_port=atoi(argv[++i]);else if(!strcmp(argv[i],"--window")&&i+1<argc)cfg.epoch_window_ms=atoi(argv[++i]);else if(!strcmp(argv[i],"--eph-host")&&i+1<argc)cfg.eph_host=argv[++i];else if(!strcmp(argv[i],"--eph-port")&&i+1<argc)cfg.eph_port=atoi(argv[++i]);else if(!strcmp(argv[i],"--ssr-host")&&i+1<argc)cfg.ssr_host=argv[++i];else if(!strcmp(argv[i],"--ssr-port")&&i+1<argc)cfg.ssr_port=atoi(argv[++i]);else if(!strcmp(argv[i],"--out-port")&&i+1<argc)cfg.out_port=atoi(argv[++i]);else if(!strcmp(argv[i],"--offset")&&i+3<argc){cfg.use_llh_offset=1;cfg.dn=atof(argv[++i]);cfg.de=atof(argv[++i]);cfg.du=atof(argv[++i]);}else if(!strcmp(argv[i],"--vbs")&&i+3<argc){cfg.use_llh_offset=0;cfg.vbs_lat=atof(argv[++i]);cfg.vbs_lon=atof(argv[++i]);cfg.vbs_alt=atof(argv[++i]);}else if(!strcmp(argv[i],"--trop"))cfg.apply_trop=1;else if(!strcmp(argv[i],"--iono"))cfg.iono_interp=1;else if(!strcmp(argv[i],"--msm")&&i+1<argc)cfg.msm_out_type=atoi(argv[++i]);else if(!strcmp(argv[i],"--no-eph"))cfg.eph_out_enable=0;else if(!strcmp(argv[i],"--no-ssr"))cfg.ssr_passthrough=0;else if(!strcmp(argv[i],"--trace")&&i+1<argc)cfg.trace_file=argv[++i];else if(!strcmp(argv[i],"--level")&&i+1<argc)cfg.trace_level=atoi(argv[++i]);else{usage(argv[0]);return 1;}}if(cfg.msm_out_type!=4&&cfg.msm_out_type!=7){fprintf(stderr,"[cfg] --msm must be 4 or 7\n");return 1;}if(cfg.iono_interp&&!cfg.dual){fprintf(stderr,"[cfg] --iono currently requires --dual\n");return 1;}setup_signals();fprintf(stderr,"rtcm_vbs_hp: A=%s:%d -> out=:%d | mode=%s | trop=%d | iono=%d | MSM=%d | SSR fwd=%d\n",cfg.src_host,cfg.src_port,cfg.out_port,cfg.use_llh_offset?"NED-offset":"absolute-LLH",cfg.apply_trop,cfg.iono_interp,cfg.msm_out_type,cfg.ssr_passthrough);run(&cfg);return 0;}
