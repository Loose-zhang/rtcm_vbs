/*------------------------------------------------------------------------------
 * tcp_io.h : minimal portable TCP client + broadcast-server for rtcm_vbs_hp
 *----------------------------------------------------------------------------*/
#ifndef TCP_IO_H
#define TCP_IO_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- TCP client (blocking, auto-reconnect in main loop) ------------------ */
typedef struct {
    int  fd;
    char host[128];
    int  port;
} tcp_client_t;

int  tcpc_connect(tcp_client_t *c, const char *host, int port);
int  tcpc_read   (tcp_client_t *c, void *buf, size_t n);  /* >0 bytes / 0 EOF / -1 err */
void tcpc_close  (tcp_client_t *c);

/* ---- TCP broadcast server (1..N simultaneous receivers) ------------------ */
typedef struct tcp_server_s tcp_server_t;

tcp_server_t *tcps_create(int port, int max_clients);
void          tcps_destroy(tcp_server_t *s);
/* non-blocking: never stalls the VBS pipeline; dropped clients are pruned. */
void          tcps_broadcast(tcp_server_t *s, const void *buf, size_t n);
int           tcps_num_clients(tcp_server_t *s);

#ifdef __cplusplus
}
#endif
#endif /* TCP_IO_H */
