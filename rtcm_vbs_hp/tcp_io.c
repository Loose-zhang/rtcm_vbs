/*------------------------------------------------------------------------------
 * tcp_io.c : POSIX (macOS/Linux) + Winsock TCP helpers
 *----------------------------------------------------------------------------*/
#include "tcp_io.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef _WIN32
  #include <winsock2.h>
  #include <ws2tcpip.h>
  #include <windows.h>
  #pragma comment(lib, "ws2_32.lib")
  typedef int socklen_t;
  #define SOCK_ERR(x) ((x) == SOCKET_ERROR)
  #define close_sock closesocket
  #define SET_NONBLOCK(fd) do { u_long m=1; ioctlsocket((fd), FIONBIO, &m); } while(0)
  #define MTX_INIT(m)    InitializeCriticalSection(m)
  #define MTX_DESTROY(m) DeleteCriticalSection(m)
  #define MTX_LOCK(m)    EnterCriticalSection(m)
  #define MTX_UNLOCK(m)  LeaveCriticalSection(m)
#else
  #include <sys/types.h>
  #include <sys/socket.h>
  #include <arpa/inet.h>
  #include <netinet/in.h>
  #include <netinet/tcp.h>
  #include <unistd.h>
  #include <fcntl.h>
  #include <errno.h>
  #define SOCK_ERR(x) ((x) < 0)
  #define close_sock close
  #define SET_NONBLOCK(fd) fcntl((fd), F_SETFL, fcntl((fd), F_GETFL, 0) | O_NONBLOCK)
  #define MTX_INIT(m)    pthread_mutex_init(m, NULL)
  #define MTX_DESTROY(m) pthread_mutex_destroy(m)
  #define MTX_LOCK(m)    pthread_mutex_lock(m)
  #define MTX_UNLOCK(m)  pthread_mutex_unlock(m)
#endif

static void tcp_init(void)
{
#ifdef _WIN32
    static int done = 0;
    if (!done) { WSADATA w; WSAStartup(MAKEWORD(2,2), &w); done = 1; }
#endif
}

/* ---- client -------------------------------------------------------------- */
int tcpc_connect(tcp_client_t *c, const char *host, int port)
{
    tcp_init();
    c->fd = -1;
    strncpy(c->host, host, sizeof(c->host)-1);
    c->host[sizeof(c->host)-1] = 0;
    c->port = port;

    int s = (int)socket(AF_INET, SOCK_STREAM, 0);
    if (s < 0) return -1;

    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port   = htons(port);
    if (inet_pton(AF_INET, host, &addr.sin_addr) != 1) {
        close_sock(s); return -1;
    }
    if (connect(s, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        close_sock(s); return -1;
    }

    int one = 1;
    setsockopt(s, IPPROTO_TCP, TCP_NODELAY, (const char*)&one, sizeof(one));
    c->fd = s;
    return 0;
}

int tcpc_read(tcp_client_t *c, void *buf, size_t n)
{
    if (c->fd < 0) return -1;
    int r = (int)recv(c->fd, (char*)buf, (int)n, 0);
    if (r == 0) return 0;
    if (r < 0) return -1;
    return r;
}

void tcpc_close(tcp_client_t *c)
{
    if (c->fd >= 0) { close_sock(c->fd); c->fd = -1; }
}

/* ---- server -------------------------------------------------------------- */
#ifndef _WIN32
#include <pthread.h>
#endif

struct tcp_server_s {
    int listen_fd;
    int port;
    int max_clients;
    int *clients;
    int  n_clients;
#ifdef _WIN32
    CRITICAL_SECTION mtx;
    HANDLE           accept_thr;
#else
    pthread_mutex_t mtx;
    pthread_t       accept_thr;
#endif
    volatile int    running;
};

#ifdef _WIN32
static DWORD WINAPI accept_thread_win(void *arg)
#else
static void *accept_thread(void *arg)
#endif
{
    tcp_server_t *s = (tcp_server_t*)arg;
    while (s->running) {
        struct sockaddr_in cli;
        socklen_t sl = sizeof(cli);
        int fd = (int)accept(s->listen_fd, (struct sockaddr*)&cli, &sl);
        if (fd < 0) {
#ifndef _WIN32
            if (errno == EINTR) continue;
#endif
            break;
        }
        int one = 1;
        setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, (const char*)&one, sizeof(one));
        SET_NONBLOCK(fd);

        MTX_LOCK(&s->mtx);
        if (s->n_clients < s->max_clients) {
            s->clients[s->n_clients++] = fd;
            char ip[64]; inet_ntop(AF_INET, &cli.sin_addr, ip, sizeof(ip));
            fprintf(stderr, "[tcp] client connected %s:%d (total=%d)\n",
                    ip, ntohs(cli.sin_port), s->n_clients);
        } else {
            close_sock(fd);
        }
        MTX_UNLOCK(&s->mtx);
    }
#ifdef _WIN32
    return 0;
#else
    return NULL;
#endif
}

tcp_server_t *tcps_create(int port, int max_clients)
{
    tcp_init();
    tcp_server_t *s = (tcp_server_t*)calloc(1, sizeof(*s));
    if (!s) return NULL;
    s->max_clients = max_clients > 0 ? max_clients : 8;
    s->clients = (int*)calloc(s->max_clients, sizeof(int));
    s->port    = port;
    MTX_INIT(&s->mtx);

    int fd = (int)socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) { free(s->clients); free(s); return NULL; }
    int one = 1;
    setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, (const char*)&one, sizeof(one));

    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family      = AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_ANY);
    addr.sin_port        = htons(port);
    if (bind(fd, (struct sockaddr*)&addr, sizeof(addr)) < 0 ||
        listen(fd, 4) < 0) {
        close_sock(fd); free(s->clients); free(s);
        return NULL;
    }
    s->listen_fd = fd;
    s->running   = 1;
#ifdef _WIN32
    s->accept_thr = CreateThread(NULL, 0, accept_thread_win, s, 0, NULL);
    if (!s->accept_thr) {
        close_sock(fd); MTX_DESTROY(&s->mtx); free(s->clients); free(s);
        return NULL;
    }
#else
    if (pthread_create(&s->accept_thr, NULL, accept_thread, s) != 0) {
        close_sock(fd); MTX_DESTROY(&s->mtx); free(s->clients); free(s);
        return NULL;
    }
#endif
    fprintf(stderr, "[tcp] VBS output listening on 0.0.0.0:%d\n", port);
    return s;
}

void tcps_destroy(tcp_server_t *s)
{
    if (!s) return;
    s->running = 0;
    close_sock(s->listen_fd);
#ifdef _WIN32
    WaitForSingleObject(s->accept_thr, INFINITE);
    CloseHandle(s->accept_thr);
#else
    pthread_join(s->accept_thr, NULL);
#endif
    MTX_LOCK(&s->mtx);
    for (int i = 0; i < s->n_clients; i++) close_sock(s->clients[i]);
    MTX_UNLOCK(&s->mtx);
    MTX_DESTROY(&s->mtx);
    free(s->clients);
    free(s);
}

void tcps_broadcast(tcp_server_t *s, const void *buf, size_t n)
{
    if (!s || n == 0) return;
    MTX_LOCK(&s->mtx);
    int w = 0;
    for (int i = 0; i < s->n_clients; i++) {
        int r = (int)send(s->clients[i], (const char*)buf, (int)n,
#ifdef MSG_NOSIGNAL
                          MSG_NOSIGNAL
#else
                          0
#endif
                          );
        if (r < 0) {
            close_sock(s->clients[i]);
            continue;   /* drop */
        }
        s->clients[w++] = s->clients[i];
    }
    s->n_clients = w;
    MTX_UNLOCK(&s->mtx);
}

int tcps_num_clients(tcp_server_t *s)
{
    if (!s) return 0;
    MTX_LOCK(&s->mtx);
    int n = s->n_clients;
    MTX_UNLOCK(&s->mtx);
    return n;
}
