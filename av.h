#ifndef AV_LIB_H
#define AV_LIB_H

#define PLATFORM_LINUX 1
#define PLATFORM_APPLE 2
#define PLATFORM_WINDOWS 3

#ifdef __linux__
    #define AV_PLATFORM PLATFORM_LINUX
#endif

#ifdef _WIN32
    #define AV_PLATFORM PLATFORM_WINDOWS
#endif

#ifdef __APPLE__
    #define AV_PLATFORM PLATFORM_APPLE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdarg.h>
#include <assert.h>

//create custom defineS for platforms?????????
#if AV_PLATFORM == PLATFORM_WINDOWS
#include <Windows.h>
#include <winsock2.h>
#pragma comment( lib, "wsock32.lib")
#else

#include <sys/stat.h>
#include <dlfcn.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <fcntl.h>
//#include <netdb.h>

#endif  //platform

typedef uint8_t         ubyte;
typedef int8_t          byte;
typedef int8_t          int8;
typedef int16_t         int16;
typedef int32_t         int32;
typedef int64_t         int64;
typedef uint8_t         uint8;
typedef uint16_t        uint16;
typedef uint32_t        uint32;
typedef uint64_t        uint64; 

typedef int8_t          i8;
typedef int16_t         i16;
typedef int32_t         i32;
typedef int64_t         i64;
typedef uint8_t         u8;
typedef uint16_t        u16;
typedef uint32_t        u32;
typedef uint64_t        u64; 

#define Kilobytes(x) ((x)*1024LL)
#define Megabytes(x) (Kilobytes(x)*1024LL)
#define Gigabytes(x) (Megabytes(x)*1024LL)

int     checkArg(int argc, char** argv, const char* argument);
int     copyFile(const char* source, const char* destination);
uint64  fileChangeTime(const char* path);
int     strToInt(const char* str);

int     f_randi(int index);

void*   loadLibrary(const char* file_name);
int     freeLibrary(void* library);
void*   getSymAddress(void* library, const char* sym_name);

void*   virtualAlloc(void* address, uint64 size);
int     virtualFree(void* address, uint64 size);

//network shit
int     initSockets();
void    quitSockets();


//IPaddress contains ip and port in network byte order | should it be in network byte order?
struct IPaddress
{
    union{
        uint32  ip;
        ubyte   octaves[4];
    };

    union{
        uint16  port;
        ubyte   port_bytes[2];
    };
};

struct Socket
{
    int _handle;
    
    int create(uint16 port)
    {
        _handle = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP); 
        if(_handle <= 0)
           return -1; 

        sockaddr_in address = {};
        address.sin_family = AF_INET;
        address.sin_addr.s_addr = INADDR_ANY;
        address.sin_port = htons(port);
        if( bind(_handle, (struct sockaddr*)&address, sizeof(address)) < 0)
            return -1;
        
        return 0;
    }

    void destroy()
    {
    #if AV_PLATFORM == PLATFORM_WINDOWS
        Closesocket(_handle); 
    #else 
        close(_handle);
    #endif
    }

    int send(void* data, int size, IPaddress ip)
    {
        sockaddr_in addr = {};
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = ip.ip;
        addr.sin_port = ip.port;
        return sendto(_handle, data, size, 0, (const sockaddr*)&addr, sizeof(addr));
    }

    int receive(void* data, int size, IPaddress* ip)
    {
        sockaddr_in addr = {};
        socklen_t len = sizeof(addr);
        int rec = recvfrom(_handle, data, size, 0, (sockaddr*)&addr, &len);
        ip->ip = addr.sin_addr.s_addr;
        ip->port = addr.sin_port;
        return rec;
    }

    int setNonBlocking(bool nonblock)
    {
    #if AV_PLATFORM == PLATFORM_WINDOWS
        int mode = nonblock;
        return ioctlsocket(_handle, FIONBIO, &mode);
    #else
        int opts = fcntl(_handle, F_GETFL);
        if(nonblock)
            opts |= O_NONBLOCK;
        else 
            opts &= ~O_NONBLOCK; 
        return fcntl(_handle, F_SETFL, opts);
    #endif
    }
};


#ifdef AV_LIB_IMPLEMENTATION

int f_randi(int index)
{
    index = (index << 13)^index;
    return ((index * (index * index * 15731 + 789221) + 1376312589) & 0x7FFFFFFF);
}

int checkArg(int argc, char** argv, const char* argument)
{
    for(int i = 0; i < argc; ++i){
        if(strcmp(argv[i], argument) == 0)
            return i;
    }
    return -1;
}

//#ifdef  PLATFORM_LINUX 

int copyFile(const char* source, const char* destination)
{
#if AV_PLATFORM == PLATFORM_WINDOWS
    int result = CopyFile(src, dst, FALSE);
    if(result == 0){
        return -1;
    }
    return 0;
#else
    //printf("copying %s to %s\n", source, destination);
    FILE* src = fopen(source, "rb");
    FILE* dst = fopen(destination, "wb+");
    int status = 0;
    uint64 size; 
    uint64 bytes_copyed = 0;
    ubyte* buffer;

    if(src == NULL || dst == NULL){
        status = -1;
        goto exit;
    }

    fseek(src, 0, SEEK_END);
    size = ftell(src); 
    fseek(src, 0, SEEK_SET);
    buffer = (ubyte*)malloc(size);
    if(buffer == NULL){
        status = -1; 
        goto exit;
    }

    fread(buffer, sizeof(ubyte), size, src);
    bytes_copyed = fwrite(buffer, sizeof(ubyte), size, dst);
    
    if(bytes_copyed != size){
        printf("wrong size \n");
        status = -1;
    }

exit:
    fclose(src);
    fclose(dst);
    if(buffer) free(buffer);
    //printf("done copying %d original size%ld size:%ld\n\n", status, size, bytes_copyed);
    return status;
#endif
}

uint64 fileChangeTime(const char* path)
{
#if AV_PLATFORM == PLATFORM_WINDOWS
    WIN32_FIND_DATA find_data;
    HANDLE file_handle = FindFirstFileA(path, &find_data);
    if(file_handle == INVALID_HANDLE_VALUE)
        return 0;
    FindClose(file_handle);
    FILETIME file_time = find_data.ftLastWriteTime;

    uint64 result = file_time.dwHighDateTime;
    result = result << 32;
    result |= file_time.dwLowDateTime;
    return result;
#else
    struct stat attr;    
    if(stat(path, &attr) == -1)
        return 0;
    return attr.st_mtime;
#endif
}

void* loadLibrary(const char* file_name)
{
#if AV_PLATFORM == PLATFORM_WINDOWS
    return LoadLibrary(file_name);
#else
    return dlopen(file_name, RTLD_NOW);
#endif
}

int freeLibrary(void* library)
{
#if AV_PLATFORM == PLATFORM_WINDOWS
    int result = FreeLibrary((HMODULE)library);
    if(result != 0) //success
        return 0;
    else 
        return -1; 
#else
    return dlclose(library);
#endif
}

void* getSymAddress(void* library, const char* sym_name)
{
#if AV_PLATFORM == PLATFORM_WINDOWS
    void* lel = GetProcAddress((HMODULE)library, sym_name);
    return lel;
#else
    return dlsym(library, sym_name);
#endif
}

void* virtualAlloc(void* address, uint64 size)
{
#if AV_PLATFORM == PLATFORM_WINDOWS
    return VirtualAlloc(address, size, MEM_COMMIT|MEM_RESERVE, PAGE_READWRITE);
#else
    //TODO: implement linux version
    assert(malloc(size));
    return 0;
#endif
}

int virtualFree(void* address, uint64 size)
{
#if AV_PLATFORM == PLATFORM_WINDOWS
    if(VirtualFree(address, size, MEM_RELEASE) == 0)
        return -1;
    else return 0;
#else
    //TODO: implement linux version
    //assert(free(size));
    free(address);
    return 0;
#endif
}

/* NETWORK STUFF */
int initSockets()
{ 
#if AV_PLATFORM == PLATFORM_WINDOWS
    WSADATA wsadata;
    if(WSAStartup(MAKEWORD(2, 2), &wsadata) == NO_ERROR)
        return 0;
    else return -1;
#else
    return 0;
#endif
}
void quitSockets()
{
#if AV_PLATFORM == PLATFORM_WINDOWS
    WSACleanup();
#else
#endif
}

int strToInt(const char* str) //unsafe? callers problem to pass \0 strs?
{
    //TODO: add base 16 support 0x...
    int result = 0;
    for(const char *c = str; *c != '\0'; ++c){
        if(*c >= '0' && *c <= '9'){
            result *= 10;
            result += *c - '0';
        } 
    }
    
    return result;
}
#endif //AV_LIB_IMPLEMENTATION
#endif //AV_LIB_H
