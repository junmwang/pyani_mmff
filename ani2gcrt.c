# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# define MAXCHAR 128
# define COLORTEXT "YES"
# define MAXATOM 512
# define debug 0

typedef struct {
	int  atomic_num;
	char element[5];
	char name[10];
	double x;
	double y;
	double z;
} ATOM;
ATOM atom[MAXATOM];
int natom = 0;
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char gfilename[MAXCHAR];
char line[MAXCHAR];
char molname[MAXCHAR]="MOL";
int i_input = 0;
int i_output= 0;
int i_gcrt = 0;
FILE *fpin, *fpout, *fpgcrt;

int max_force=0;
int max_displacement=0;
int rms_force=0;
int rms_displacement=0;
double energy = 0;

void read() {
int i;
char tmpc1[MAXCHAR];
char tmpc2[MAXCHAR];
char tmpc3[MAXCHAR];
char tmpc4[MAXCHAR];
char tmpc5[MAXCHAR];
char tmpc6[MAXCHAR];
char tmpc7[MAXCHAR];
char tmpc8[MAXCHAR];
char tmpc9[MAXCHAR];
char tmpc10[MAXCHAR];
char tmpc11[MAXCHAR];
for (;;) {
	if (fgets(line, MAXCHAR, fpin) == NULL) {
        	break;
	}
	strcpy(tmpc1, "");
	strcpy(tmpc2, "");
	strcpy(tmpc3, "");
	strcpy(tmpc4, "");
	strcpy(tmpc5, "");
	strcpy(tmpc6, "");
	strcpy(tmpc7, "");
	strcpy(tmpc8, "");
	strcpy(tmpc9, "");
	strcpy(tmpc10, "");
	strcpy(tmpc11, "");

	tmpc1[0]='\0';
	tmpc2[0]='\0';
	tmpc3[0]='\0';
	tmpc4[0]='\0';
	tmpc5[0]='\0';
	tmpc6[0]='\0';
	tmpc7[0]='\0';
	tmpc8[0]='\0';
	tmpc9[0]='\0';
	tmpc10[0]='\0';
	tmpc11[0]='\0';
	sscanf(line, "%s%s%s%s%s%s%s%s%s%s%s", tmpc1, tmpc2, tmpc3, tmpc4, tmpc5, tmpc6, tmpc7, tmpc8, tmpc9, tmpc10, tmpc11);
	if(strcmp(tmpc1, "Maximum") == 0 && strcmp(tmpc2, "Force") == 0) {
		if(strcmp(tmpc5, "YES") == 0) max_force = 1;
	}
	if(strcmp(tmpc1, "RMS") == 0 && strcmp(tmpc2, "Force") == 0) {
		if(strcmp(tmpc5, "YES") == 0) rms_force = 1;
	}
	if(strcmp(tmpc1, "Maximum") == 0 && strcmp(tmpc2, "Displacement") == 0) {
		if(strcmp(tmpc5, "YES") == 0) max_displacement = 1;
	}
	if(strcmp(tmpc1, "RMS") == 0 && strcmp(tmpc2, "Displacement") == 0) {
		if(strcmp(tmpc5, "YES") == 0) rms_displacement = 1;
	}
	if(strcmp(tmpc1, "Energy:") == 0) energy=atof(tmpc2);
	
	if(strcmp(tmpc1, "ATOM") == 0) {
		strcpy(atom[natom].name, tmpc3);
		atom[natom].x= atof(tmpc6);
		atom[natom].y= atof(tmpc7);
		atom[natom].z= atof(tmpc8);
		strcpy(atom[natom].element, tmpc11);
		natom++;
	}
 }
	if(debug) {
		for(i=0; i<natom; i++) 
			printf("%4d %5s %5s %9.4lf %9.4lf %9.4lf\n", i+1, atom[i].name, atom[i].element, atom[i].x, atom[i].y,atom[i].z);
	}
}

void write1() {
int i;
fprintf(fpout, "--Link1--\n");
fprintf(fpout, "%s\n", "%nproc=2");
fprintf(fpout, "%s\n", "%mem=2GB");
fprintf(fpout, "%s%s\n", "%chk=",molname);
fprintf(fpout, "%s\n", "#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt");
fprintf(fpout, "\nRemark line goes here\n\n0   1\n"); 

for(i=0; i<natom;i++) {
	fprintf(fpout, "%5s %18.10lf %18.10lf %18.10lf\n", atom[i].element, atom[i].x, atom[i].y, atom[i].z);
}
fprintf(fpout, "\n\n");
}

void write2() {
int i;
int space=0;
int flag=1;

for (;;) {
	if (fgets(line, MAXCHAR, fpgcrt) == NULL) { break; }
	if(strlen(line) <= 2) {
		space++;
		fprintf(fpout, line);
		continue;
	}
	if(space == 2 && flag == 1) {
		fprintf(fpout, "%s", line);
		for(i=0; i<natom;i++) {
			fprintf(fpout, "%5s %18.10lf %18.10lf %18.10lf\n", 
				atom[i].element, atom[i].x, atom[i].y, atom[i].z);
		}	
		flag = 0;
		continue;
	}	
	if(space == 2 && flag == 0) continue;
	if(space < 2 || space >= 3) fprintf(fpout, "%s", line);
}

}

int main(int argc, char *argv[]) {
	int i, flag;
	int count = 0;
	int sid=0;
	char cmdstr[MAXCHAR];
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("[31mUsage: ani2gcrt -i [0m ani logfile \n"
				 "[31m                -o [0m output file \n"
				 "[31m                -g [0m gcrt input file, optional \n");
			exit(0);
		}
		if (argc != 5 && argc != 7) {
			printf
				("[31mUsage: ani2gcrt -i [0m ani logfile \n"
				 "[31m                -o [0m output file \n"
				 "[31m                -g [0m gcrt input file, optional \n");
			exit(0);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("Usage: ani2gcrt -i  ani logfile \n"
				 "                -o  output file \n"
				 "                -g  gcrt input file, optional \n");
			exit(0);
		}
		if (argc != 5 && argc != 7) {
			printf
				("Usage: ani2gcrt -i  ani logfile \n"
				 "                -o  output file \n"
				 "                -g  gcrt input file, optional \n");
		}
	}

	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			i_input = 1;
		}
		if (strcmp(argv[i], "-o") == 0)  {
			strcpy(ofilename, argv[i + 1]);
			i_output = 1;
		}
		if (strcmp(argv[i], "-g") == 0) {
			strcpy(gfilename, argv[i + 1]);
			i_gcrt = 1;
		}
	}
	if ((fpin = fopen(ifilename, "r")) == NULL) {
        	fprintf(stderr, "Cannot open ani input file %s to read, exit\n", ifilename);
        	exit(1);
	}
	if(i_gcrt == 1) {
		if ((fpgcrt = fopen(gfilename, "r")) == NULL) {
        		fprintf(stderr, "Cannot open ani log file %s to read, exit\n", gfilename);
        		exit(1);
		}
	}
	if ((fpout = fopen(ofilename, "w")) == NULL) {
        	fprintf(stderr, "Cannot open output file %s to write, exit\n", ofilename);
        	exit(1);
	}
	read();
	count=0;
	sid = 0;
	for(i=0; i<strlen(ifilename); i++) {
		if(ifilename[i]=='/') sid = i+1;
	}
	
	for(i=sid; i<strlen(ifilename); i++) {
		if(ifilename[i] == '.') {
			molname[count]='\0';
			break;	
		}
		molname[count++]=ifilename[i];		
	}
	if(i_gcrt == 0) write1();
	if(i_gcrt == 1) write2();
fclose(fpin);
if(i_gcrt == 1) fclose(fpgcrt);
fclose(fpout);
return 0;
}
