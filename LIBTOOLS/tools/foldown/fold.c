/*
 * foldonw
 * -------
 *
 */

#include <stdio.h>

#define MAX_LINE_LENGTH 60

void foldonw ( FILE *fp, int csp, char *ssp, int pos )
{
    int c;
    int l = 0;
    int split = 0;

    while ( (c=getc(fp)) != EOF )
    {
        if ( c == '\n' )
            l = split = 0;
        else
        {
            l++;
            if ( l > pos )
                split = 1;
            if ( split && c == csp )
                split = 2;
        }
        putchar(c);
        if ( split == 2 )
        {
            printf("\n%s", ssp);
            l = split = 0;
        }
    }

    return;
}

int main ( int argc, char **argv )
{
    int   iarg = 0;
    char *sf = NULL;
    char *ssp = "";
    FILE *fp = stdin;
    int   csp = ',';
    int   pos = MAX_LINE_LENGTH;

    while ( ++iarg < argc && *(argv[iarg]++) == '-' )
        switch ( *argv[iarg] )
        {
            case 'f' : sf = argv[++iarg]; break;
            case 'p' : pos = atoi(argv[++iarg]); break;
            case 'c' : csp = *argv[++iarg]; break;
            case 's' : ssp = argv[++iarg]; break;
            default :
                fprintf(stderr, "Usage: foldonw [-f filename] [-p pos] [-c char] [-s begin-string]\n");
                exit(1);
        }

    if ( sf != NULL && (fp=fopen(sf,"r")) == NULL )
    {
        fprintf(stderr, "%s: no such file or directory\n", sf);
        exit(1);
    }

    foldonw(fp, csp, ssp, pos);

    return 0;
}
