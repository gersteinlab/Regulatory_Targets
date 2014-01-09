/*
 * vmstat.cpp
 *
 *  Created on: Sep 22, 2011
 *      Author: lousk
 */

#include "vmstat.h"

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <ctype.h>


void getVM(){

        pid_t pid;
        if((pid=getpid())<0){
        printf("error,could not get pid");
        }
        char vm[11][10]={"VmPeak","VmSize","VmLck",
        		 	   "VmHWM","VmRSS","VmData",
        		 	   "VmStk","VmExe","VmLib",
        		 	   "VmPTE","VmSwap"};

        char c;
        char pidbuf[128];
        printf("DEBUG:pid=%d ",(int)pid);
        sprintf(pidbuf,"/proc/%d/status",(int)pid);
        FILE *stat;
        stat=fopen(pidbuf,"r");
        int i=0,j,k;
        char line[1024];
        char head[1024];
        char cont[1024];
        c=fgetc(stat);
        while(!feof(stat)){
                if(c!='\n'){
                        line[i]=c;
                        i++;
                        }
                        else{
                             line[i]=0x00;
                            // printf("length:%d=>%s\n",strlen(line),line);
                             for(j=0;j<i;j++){

                            	if(line[j]==':'){
                            		strncpy(head,line,j);
                            		head[j]=0x00;
                            		if(strcmp(head,vm[0])==0 || strcmp(head,vm[1])==0 ||
                            		   strcmp(head,vm[2])==0 ||strcmp(head,vm[3])==0 ||
                            		   strcmp(head,vm[4])==0 ||strcmp(head,vm[5])==0 ||
                            		   strcmp(head,vm[6])==0 ||strcmp(head,vm[11])==0
                            		){
                            		j++;
                            		while(isspace(line[j])){
                            			j++;
                            			}
                            		strncpy(cont,line+j,(i-j));
                            		cont[i-j]=0x00;
                            		printf("%s=>%s\t",head,cont);
                            		}
                            	}
                             }

                             //printf("\n>>>>%s=>%s\n",head,cont);
                             i=0;
                         }

                c=fgetc(stat);
        }

        //delete[] pidbuf;
        //delete[] line;
        //delete[] head;
        //delete[] cont;
        //delete[] vm;
        printf("\n");
}



void getVM(const char* info){

        pid_t pid;
        if((pid=getpid())<0){
        printf("error,could not get pid");
        }
        char vm[11][10]={"VmPeak","VmSize","VmLck",
        		 	   "VmHWM","VmRSS","VmData",
        		 	   "VmStk","VmExe","VmLib",
        		 	   "VmPTE","VmSwap"};

        char c;
        char pidbuf[128];
        printf("DEBUG:%s(pid:%d) ",info,(int)pid);
        sprintf(pidbuf,"/proc/%d/status",(int)pid);
        FILE *stat;
        stat=fopen(pidbuf,"r");
        int i=0,j,k;
        char line[1024];
        char head[1024];
        char cont[1024];
        c=fgetc(stat);
        while(!feof(stat)){
                if(c!='\n'){
                        line[i]=c;
                        i++;
                        }
                        else{
                             line[i]='\0';
                            // printf("length:%d=>%s\n",strlen(line),line);
                             for(j=0;j<i;j++){

                            	if(line[j]==':'){
                            		strncpy(head,line,j);
                            		head[j]=0x00;
                            		if(strcmp(head,vm[0])==0 || strcmp(head,vm[1])==0 ||
                            		   strcmp(head,vm[2])==0 ||strcmp(head,vm[3])==0 ||
                            		   strcmp(head,vm[4])==0 ||strcmp(head,vm[5])==0 ||
                            		   strcmp(head,vm[6])==0 ||strcmp(head,vm[11])==0
                            		){
                            		j++;
                            		while(isspace(line[j])){
                            			j++;
                            			}
                            		strncpy(cont,line+j,(i-j));
                            		cont[i-j]=0x00;
                            		printf("%s=>%s\t",head,cont);
                            		}
                            	}
                             }

                             //printf("\n>>>>%s=>%s\n",head,cont);
                             i=0;
                         }

                c=fgetc(stat);
        }
        printf("\n");
}


/*
int main(){

	while(1){
	getVM();
	sleep(1);
	}

	return 0;
}

*/
