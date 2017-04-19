#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#define G 0.0000005
#define MIN_DIST 1.0
const int unitSize = 8;//unitSize:每个计算单位的大小，多少个float
double wide = 1000;
double mass = 10;
double velo = 0.005;

typedef struct Orb {
  double x,y,z,vx,vy,vz,m,st;
}Orb;
typedef union {
  Orb orbst;
  double fmem[8];
}UOrb;

void initList(Orb *olist, int nUnit, int style);
void calcOne(Orb*o, int oId, Orb*olist, int nUnit);
void calcGravity(Orb*o, Orb*ta, double dist, double*gx, double*gy, double*gz);
void printList(double *list, int len, int rank, const char* s);
void saveList(Orb *olist, int len, const char* filepath);
void loadList(double **list, int *len, const char* filepath);

int main(int argc, char *argv[])
{
    int rank, value=0, nWorkers;
    int nUnit = 5, nTimes = 3;
    int ch;
    while ((ch = getopt(argc, argv, "l:m:v:w:t:")) != -1) {
        switch (ch) {
            case 'l':
                nUnit = atoi(optarg);
                printf("set nUnit=%d\n", nUnit);
                break;
            case 'm':
                mass = atof(optarg);
                printf("set mass=%f\n", mass);
                break;
            case 'v':
                velo = atof(optarg);
                printf("set velo=%f\n", velo);
                break;
            case 'w':
                wide = atof(optarg);
                printf("set wide=%f\n", wide);
                break;
            case 't':
                nTimes = atoi(optarg);
                printf("set times=%d\n", nTimes);
                break;
            defualt:
                break;
        }
    }
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nWorkers);
    // 任务如何均匀下发？
    int nUnitPerWorker = nUnit / nWorkers;//nUnit:list有多少个unit。list字节长度=sizeof(double)*unitSize
    int listSizeOfByte = sizeof(double) * unitSize * nUnit;
    int chunkSize = unitSize*nUnitPerWorker;
    double* list = (double*)malloc(listSizeOfByte);
    long calcTimes = (long)(nUnit * nUnit) * (long)nTimes;
    printf("list addr=%p nUnit=%d nUnitPerWorker=%d chunkSize=%d listbytes=%d\n", list, nUnit, nUnitPerWorker, chunkSize, listSizeOfByte);

    // init list here
    if (rank == 0) {
        initList((Orb*)list, nUnit, 0);
    }
    // 广播必须每个rank都调用，否则报 Fatal error in PMPI_Barrier: Message truncated, error stack
    MPI_Bcast(list, nUnit*unitSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    printList(list, nUnit*unitSize, rank, "after init bcast");
    double startTime = MPI_Wtime(), endTime;

    int i = 0, j = 0, k = 0;
    for (j=0; j<nTimes; ++j) {
        // calculate list here
        //list[rank*chunkSize] = (double)rank*100+(double)j;
        for (k=0; k<nUnitPerWorker; ++k) {
            Orb* o = (Orb*)(list+rank*nUnitPerWorker*unitSize+k);
            int oId = rank*nUnitPerWorker+k;
            calcOne(o, oId, (Orb*)list, nUnit);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        // bcast 对于发送方来说等于send，对于接受方来说等于recv，调用时需要for i in nWorkers来依次bcast。相当于调用nWorkers*nWorkers次。
        // 不加for循环，root设置为自己，相当于只有root.send,没有other.recv。
        for (i=0; i<nWorkers; ++i) MPI_Bcast(list+i*chunkSize, chunkSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printList(list, nUnit*unitSize, rank, "finally");
    printf("list+i*chunkSize=%p\n", list+i*chunkSize);

    endTime = MPI_Wtime();
    printf("all used time:%f calc times:%ld cps:%eg\n", endTime-startTime, calcTimes, (double)calcTimes/(endTime-startTime));
    saveList((Orb*)list, nUnit, "./thelist1");
    free(list);
    list = NULL;
    MPI_Finalize();
    return 0;
}

void printList(double *list, int len, int rank, const char* s) {
    //char str[1000] = {};
    int bufferlen = 100000;
    static char* str = NULL;
    if (str == NULL) {
        str = (char*) malloc(bufferlen);//
    } else {
        memset(str, 0, bufferlen);
    }
    int i = 0, pos = 0;
    for (i=0; i<len; i++) {
        char stmp[16] = {};
        sprintf(stmp, "%f,", list[i]);
        if (((i+1)%unitSize) == 0) {
            sprintf(stmp+strlen(stmp), "\n");
        }
        memcpy((str+pos), stmp, strlen((const char*)stmp));
        pos += strlen(stmp);
        if (pos>=bufferlen) {
            printf("str not enough! i=%d len=%d rank=%d\n", i, len, rank);
            break;
        }
    }
    printf("%s, rank=%d list={\n%s}\n", s, rank, str);
}
void saveList(Orb *olist, int len, const char* filepath) {
}
/**/
void loadList(double **list, int *len, const char* filepath) {
}
void initList(Orb *olist, int nUnit, int style) {
    int i = 0;
    srandom(time(NULL));
    for (i=0; i<nUnit; ++i) {
        //Orb* o = (Orb*)(list+i*unitSize);
        Orb* o = olist+i;
        o->x  = (double)random()/(double)RAND_MAX*wide;
        o->y  = (double)random()/(double)RAND_MAX*wide;
        o->z  = (double)random()/(double)RAND_MAX*wide;
        o->vx = (double)random()/(double)RAND_MAX*velo;
        o->vy = (double)random()/(double)RAND_MAX*velo;
        o->vz = (double)random()/(double)RAND_MAX*velo;
        o->m  = (double)random()/(double)RAND_MAX*mass;
    }
}
void calcGravity(Orb*o, Orb*ta, double dist, double*gx, double*gy, double*gz) {
    double a = ta->m / (dist*dist) * G;
    *gx = - a * (o->x - ta->x) / dist;
    *gy = - a * (o->y - ta->y) / dist;
    *gz = - a * (o->z - ta->z) / dist;
}
/* Orb update once with list */
void calcOne(Orb*o, int oId, Orb*olist, int nUnit) {
    //Orb* mest = (Orb*)*o;
    if (o->st == 0) {
        int i = 0;
        double g, gx, gy, gz, gax, gay, gaz, dist = 0;
        g = gx = gy = gz = gax = gay = gaz = 0;
        for (i=0; i<nUnit; ++i) {
            Orb* ta = (Orb*)(olist+i);
            if (oId == i || ta->st == 1) {
                continue;
            }
            dist = sqrt((double)((ta->x-o->x)*(ta->x-o->x)+(ta->y-o->y)*(ta->y-o->y)+(ta->x-o->z)*(ta->z-o->z)));
            int isTooRappid = dist*dist<(o->vx*o->vx+o->vy*o->vy+o->vz*o->vz)*10;
            if (dist<MIN_DIST || isTooRappid) {
                // crash
                o->st = 1;
                printf("one crash: oid=%d, tid=%d dist=%f isTooRappid=%d\n", oId, i, dist, isTooRappid);
                continue;
            }
            calcGravity(o, ta, dist, &gx, &gy, &gz);
            gax += gx;
            gay += gy;
            gaz += gz;
        }
//        gax = gax/nUnit;
//        gay = gay/nUnit;
//        gaz = gaz/nUnit;
        o->x += o->vx;
        o->y += o->vy;
        o->z += o->vz;
        o->vx += gax;
        o->vy += gay;
        o->vz += gaz;
    }
}


