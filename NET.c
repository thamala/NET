/*
 
 Network Enrichment Test (NET)
 
 Version 2019.05.19
 
 Copyright (C) 2019 Tuomas Hamala
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 For any other enquiries send an email to tuomas.hamala@gmail.com
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#define merror "\nERROR: System out of memory\n\n"
#define version "2019.05.19"

typedef struct{
    int obs, est_i;
    double exp, est, pval;
    char *name, **genes;
}Module_s;

typedef struct{
    int found;
    double est;
    char *name;
}Gene_s;

void openFiles(int argc, char *argv[]);
Gene_s *readGenes(FILE *gene_file, int *n, int t);
Module_s *readModules(FILE *module_file, int *m, int *g);
double *readNeutral(FILE *neutral_file, int *n);
void testCand(Module_s *modules, Gene_s *genes, int gene_n, int module_n, int mgene_n);
void testAll(Module_s *modules, Gene_s *genes, double neutral[], int gene_n, int module_n, int mgene_n, int neut_n, int side, double perm);
double hgProb(double facs[], int a, int b, int c, int d);
double estPval(int df, double chi);
static double igf(double s, double z);
int isNumeric (const char *s);
void lineTerminator(char *line);
void printInfo(void);

int main(int argc, char *argv[]){
    
    int second=0, minute=0, hour=0;
    time_t timer=0;
    
    if(argc == 1){
        printInfo();
        exit(EXIT_FAILURE);
    }
    
    timer = time(NULL);
    
    openFiles(argc, argv);
    
    second = time(NULL) - timer;
    
    minute = second / 60;
    
    hour = second / 3600;
    
    if(hour > 0)
        fprintf(stderr,"\nElapsed time: %i h, %i min & %i sec\n\n", hour, minute-hour*60, second-minute*60);
    else if(minute > 0)
        fprintf(stderr,"\nElapset time: %i min & %i sec\n\n", minute, second-minute*60);
    else if(second > 5)
        fprintf(stderr,"\nElapsed time: %i sec\n\n", second);
    else
        fprintf(stderr,"\n");
    
    return 0;
}

void openFiles(int argc, char *argv[]){
    
    int i, test=-1, side=0, gene_n=0, module_n=0, mgene_n=0, neut_n=0;
    double *neutral, perm=0;
    FILE *module_file=NULL, *gene_file=NULL, *neutral_file=NULL;
    Gene_s *genes;
    Module_s *modules;
    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--modules") == 0){
            if((module_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr, "\nERROR: Cannot open file %s\n", argv[i]);
                printInfo();
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t--modules %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--in") == 0){
            if((gene_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr, "\nERROR: Cannot open file %s\n", argv[i]);
                printInfo();
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t--in %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--neutral") == 0){
            if((neutral_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr, "\nERROR: Cannot open file %s\n", argv[i]);
                printInfo();
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t--neutral %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--test") == 0){
            if(isNumeric(argv[++i]))
                test = atoi(argv[i]);
            else
                test = 2;
            if(test != 0 && test != 1){
                fprintf(stderr,"\nERROR: Value for -t must be either 0 or 1\n");
                printInfo();
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t--test %i\n", test);
        }
        
        else if(strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--side") == 0){
            if(isNumeric(argv[++i]))
                side = atoi(argv[i]);
            else
                side = 2;
            if(side != 0 && side != 1){
                fprintf(stderr,"\nERROR: Value for -s must be either 0 or 1\n");
                printInfo();
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t--side %i\n", side);
        }
        
        else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--cycles") == 0){
            if(isNumeric(argv[++i]))
                perm = atof(argv[i]);
            else{
                fprintf(stderr,"\nERROR: Value for -c must be an integer\n");
                printInfo();
                exit(EXIT_FAILURE);
            }
            if(perm >= 1e5)
                fprintf(stderr, "\t--cycles %.0e\n", perm);
            else
                fprintf(stderr, "\t--cycles %.0f\n", perm);
        }
        
        else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0){
            fprintf(stderr, "\t--help\n");
            printInfo();
            exit(EXIT_FAILURE);
        }
        
        else{
            fprintf(stderr,"\nERROR: Unknown argument '%s'\n", argv[i]);
            printInfo();
            exit(EXIT_FAILURE);
        }
    }
    
    fprintf(stderr,"\n");
    
    if(module_file == NULL || gene_file == NULL){
        fprintf(stderr,"ERROR: Two input files are required\n");
        printInfo();
        exit(EXIT_FAILURE);
    }
    
    if(test == -1){
        fprintf(stderr,"ERROR: The choise of test, -t [0,1], is required\n");
        printInfo();
        exit(EXIT_FAILURE);
    }
    
    if(test == 0 && perm > 0)
        fprintf(stderr,"Warning: Permutations are ignored when performing a candidate gene test\n\n");
    
    if(test == 0 && side == 1)
        fprintf(stderr,"Warning: -s is ignored when performing a candidate gene test\n\n");
    
    if(neutral_file != NULL && perm == 0){
        fprintf(stderr,"Warning: Neutral data provided without defining permutation cycles\nUsing -c 10000\n\n");
        perm = 10000;
    }
    
    if(test == 0)
        genes = readGenes(gene_file, &gene_n, 0);
    else
        genes = readGenes(gene_file, &gene_n, 1);
    
    modules = readModules(module_file, &module_n, &mgene_n);
    
    if(neutral_file != NULL)
        neutral = readNeutral(neutral_file, &neut_n);
    
    if(test == 0)
        testCand(modules, genes, gene_n, module_n, mgene_n);
    else
        testAll(modules, genes, neutral, gene_n, module_n, mgene_n, neut_n, side, perm);
}

Gene_s *readGenes(FILE *gene_file, int *n, int t){
    
    int i, char_i=0, maxchar=0, line_i=0;
    char c, *line, *temp=NULL;
    Gene_s *list;
    
    while((c=fgetc(gene_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(gene_file);
    
    if((line = malloc((maxchar+2)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    if((list = malloc(line_i*sizeof(Gene_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<line_i;i++){
        if((list[i].name = malloc(maxchar*sizeof(char))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    i = 0;
    
    while(fgets(line, maxchar+2, gene_file) != NULL){
        if(line[0] != '#'){
            lineTerminator(line);
            temp = strtok(line, "\t");
            if(temp == NULL){
                fprintf(stderr,"ERROR: Gene-file is not in the specified format\n");
                printInfo();
                exit(EXIT_FAILURE);
            }
            else{
                strcpy(list[i].name, temp);
                list[i].found = 0;
            }
            temp = strtok(NULL, "\t");
            if(t == 1){
                if(temp == NULL){
                    fprintf(stderr,"ERROR: Gene-file must contain test statistics when -t is 1\n");
                    printInfo();
                    exit(EXIT_FAILURE);
                }
                else if(isNumeric(temp) == 0){
                    fprintf(stderr,"ERROR: Test statistic must be numeric\n");
                    printInfo();
                    exit(EXIT_FAILURE);
                }
                else
                    list[i].est = atof(temp);
            }
            else{
                if(temp != NULL && i == 0)
                    fprintf(stderr,"Warning: Gene-file contains more than one column. Did you mean to use -t 1?\n\n");
            }
            i++;
        }
    }
    
    *n = i;
    
    return list;
}

Module_s *readModules(FILE *module_file, int *m, int *g){
    
    int i, j, k, char_i=0, line_i=0, maxchar=0, mod_i=0;
    char c, **genes, **modules, *line, *temp=NULL;
    Module_s *list;
    
    while((c=fgetc(module_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(module_file);
    
    if((line = malloc((maxchar+2)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    if((genes = malloc((line_i+2)*sizeof(char*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    if((modules = malloc((line_i+2)*sizeof(char*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<line_i+2;i++){
        
        if((genes[i] = malloc(maxchar*sizeof(char))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
        
        if((modules[i] = malloc(maxchar*sizeof(char))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    line_i = 0;
    
    while(fgets(line, maxchar+2, module_file) != NULL){
        if(line[0] != '#'){
            lineTerminator(line);
            temp = strtok(line, "\t");
            if(temp == NULL){
                fprintf(stderr,"ERROR: Module-file is not in the specified format\n");
                printInfo();
                exit(EXIT_FAILURE);
            }
            else
                strcpy(genes[line_i], temp);
            temp = strtok(NULL, "\t");
            if(temp == NULL){
                fprintf(stderr,"ERROR: Module-file is not in the specified format\n");
                printInfo();
                exit(EXIT_FAILURE);
            }
            else
                strcpy(modules[line_i], temp);
            line_i++;
        }
    }
    
    for(i=0;i<line_i;i++){
        for(j=0;j<i;j++){
            if(strcmp(modules[i], modules[j]) == 0)
                break;
        }
        if(j == i)
            mod_i++;
    }
    
    if((list = malloc(mod_i*sizeof(Module_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<mod_i;i++){
        
        if((list[i].name = malloc(maxchar*sizeof(char))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
        
        if((list[i].genes = malloc(line_i*sizeof(char*))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
        
        for(j=0;j<line_i;j++){
            if((list[i].genes[j] = malloc(maxchar*sizeof(char))) == NULL){
                fprintf(stderr,merror);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    mod_i = 0;
    
    for(i=0;i<line_i;i++){
        for(j=0;j<i;j++){
            if(strcmp(modules[i], modules[j]) == 0)
                break;
        }
        if(j == i){
            strcpy(list[mod_i].name, modules[i]);
            list[mod_i].obs = 0;
            list[mod_i].est = 0;
            list[mod_i].pval = 0;
            list[mod_i].est_i = 0;
            mod_i++;
        }
    }
    
    for(i=0;i<line_i;i++){
        for(j=0;j<mod_i;j++){
            if(strcmp(modules[i], list[j].name) == 0){
                while(list[j].genes[k][0]!='\0')
                    k++;
                strcpy(list[j].genes[k], genes[i]);
                k = 0;
            }
        }
    }
    
    *m = mod_i;
    *g = line_i;
    
    for(i=0;i<line_i+2;i++){
        free(genes[i]);
        free(modules[i]);
    }
    
    free(genes);
    free(modules);
    free(line);
    fclose(module_file);
    
    return list;
}

double *readNeutral(FILE *neutral_file, int *n){
    
    int i, char_i=0, line_i=0, maxchar=0;
    double *list;
    char c, *line;
    
    while((c=fgetc(neutral_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(neutral_file);
    
    if((line = malloc((maxchar+2)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    if((list = malloc(line_i*sizeof(double))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    i = 0;
    
    while(fgets(line, maxchar+2, neutral_file) != NULL){
        lineTerminator(line);
        if(isNumeric(line) == 0){
            fprintf(stderr,"ERROR: Neutral estimates must be numeric\n");
            printInfo();
            exit(EXIT_FAILURE);
        }
        else
            list[i] = atof(line);
        i++;
    }
    
    *n = i;
    
    return list;
}

void testCand(Module_s *modules, Gene_s *genes, int gene_n, int module_n, int mgene_n){
    
    int i, j, k, found_i=0, gene_i=0, print_i=0, df=0;
    int a=0, b=0, c=0, d=0, n=0, m=0;
    double mltp=0, x=0, chi=0, pval=0, *facs;
    
    //Counts the number of candidate genes found in each module
    for(i=0;i<gene_n;i++){
        for(j=0;j<module_n;j++){
            for(k=0;modules[j].genes[k][0]!=0;k++){
                if(strcmp(genes[i].name, modules[j].genes[k]) == 0){
                    modules[j].obs++;
                    if(genes[i].found == 0){
                        print_i++;
                        genes[i].found = 1;
                    }
                    found_i++;
                }
            }
        }
    }
    
    if((facs = malloc(mgene_n*sizeof(double))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    facs[0] = 0;
    
    //Calculates log-transformed factorials for Fisher's test
    for(i=1;i<=mgene_n;i++)
        facs[i] = facs[i-1] + log((double)i);
    
    printf("Module\tObserved\tExpected\tFold\tP.value\n");
    
    for(i=0;i<module_n;i++){
        
        gene_i = 0;
        while(modules[i].genes[gene_i][0]!='\0')
            gene_i++;
        mltp = (double)gene_i / (double)mgene_n;
        modules[i].exp = (double)found_i * mltp;
        
        //Module specific p-values are estimated with one-sided Fisher's exact test
        a = modules[i].obs;
        b = gene_i - a;
        c = found_i - a;
        d = mgene_n - found_i - b;
        n = a + b + c + d;
        
        pval = 0;
        
        for(j=a;j<=n;j++){
            if(a+b-j >= 0 && a+c-j >= 0 && d-a+j >= 0)
                pval += exp(hgProb(facs,j,a+b-j,a+c-j,d-a+j));
        }
        
        printf("%s\t%i\t%.2f\t%.3f\t%f\n", modules[i].name, modules[i].obs, modules[i].exp, modules[i].obs/modules[i].exp, pval);
    }
    
    if(module_n < 100){
        //If module n < 100, estimates overall p-value with chi-square test
        for(i=0;i<module_n;i++){
            x = (double)modules[i].obs - modules[i].exp;
            chi += x * x / modules[i].exp;
        }
        
        df = module_n - 1;
        
        pval = estPval(df, chi);
    }
    
    for(i=0;i<module_n;i++){
        for(j=0;j<mgene_n;j++)
            free(modules[i].genes[j]);
        free(modules[i].genes);
        free(modules[i].name);
    }
    free(modules);
    for(i=0;i<gene_n;i++)
        free(genes[i].name);
    free(genes);
    free(facs);
    
    if(isatty(1))
        fprintf(stderr,"\n");
    
    fprintf(stderr,"Done!\n\n%i modules contain a total of %i genes\n%i out of %i genes found from modules\n", module_n, mgene_n, print_i, gene_n);
    if(module_n < 100)
        fprintf(stderr,"\nChi-square = %.2f\nDf = %i\nP-value = %f\n", chi, df, pval);
}

void testAll(Module_s *modules, Gene_s *genes, double neutral[], int gene_n, int module_n, int mgene_n, int neut_n, int side, double perm){
    
    int i, j, k, r_gene=0, found=0, found_i=0, print_i=0;
    double r_est=0, mean=0;
    
    //Calculates the summed test statistic for each module
    for(i=0;i<gene_n;i++){
        found = 0;
        for(j=0;j<module_n;j++){
            for(k=0;modules[j].genes[k][0]!=0;k++){
                if(strcmp(genes[i].name, modules[j].genes[k]) == 0){
                    modules[j].est += genes[i].est;
                    modules[j].est_i++;
                    mean += genes[i].est;
                    found_i++;
                    if(genes[i].found == 0){
                        print_i++;
                        genes[i].found = 1;
                    }
                }
            }
        }
    }
    
    //P-values are estimated with permutation test
    if(perm > 0){
        //Sets seed for permutation testing
        srand(time(NULL));
        for(i=0;i<module_n;i++){
            if(modules[i].est_i > 0){
                for(j=0;j<perm;j++){
                    for(k=0;k<modules[i].est_i;k++){
                        //Using neutral data if provided
                        if(neut_n > 0){
                            r_gene = rand() % neut_n;
                            r_est += neutral[r_gene];
                        }
                        else{
                            r_gene = rand() % gene_n;
                            if(genes[r_gene].found == 1)
                                r_est += genes[r_gene].est;
                            else
                                k--;
                        }
                    }
                    if(side == 1){
                        if(r_est <= modules[i].est)
                            modules[i].pval++;
                    }
                    else{
                        if(r_est >= modules[i].est)
                            modules[i].pval++;
                    }
                    
                    r_est = 0;
                }
                modules[i].pval /= perm;
            }
        }
    }
    
    if(perm > 0)
        printf("Module\tGenes\tMean\tP.value\n");
    else
        printf("Module\tGenes\tMean\n");
    
    for(i=0;i<module_n;i++){
        if(modules[i].est_i > 0){
            printf("%s\t%i\t%.3f\t", modules[i].name, modules[i].est_i, modules[i].est/(double)modules[i].est_i);
            if(perm > 0){
                printf("%f\n", modules[i].pval);
            }
            else
                printf("\n");
        
        }
    }
    
    for(i=0;i<module_n;i++){
        for(j=0;j<mgene_n;j++)
            free(modules[i].genes[j]);
        free(modules[i].genes);
        free(modules[i].name);
    }
    free(modules);
    for(i=0;i<gene_n;i++)
        free(genes[i].name);
    free(genes);
    
    if(isatty(1))
        fprintf(stderr,"\n");
    
    fprintf(stderr,"Done!\n\n%i modules contain a total of %i genes\n%i out of %i genes found from modules\n", module_n, mgene_n, print_i, gene_n);
    fprintf(stderr,"Mean test statistic = %.3f\n", mean/(double)found_i);
}

double hgProb(double facs[], int a, int b, int c, int d){
    
    return facs[a+b] + facs[c+d] + facs[a+c] + facs[b+d] - facs[a] - facs[b] - facs[c] - facs[d] - facs[a+b+c+d];
}

double estPval(int df, double chi){
    
    double k=0, x=0, p=0;
    
    if(chi < 0 || df < 1)
        return 1;
    
    k = (double)df * 0.5;
    x = chi * 0.5;
    
    if(df == 2)
        return exp(-1 * x);
    
    p = igf(k, x);
    
    if(isnan(p) || isinf(p) || p <= 1e-8)
        return 1;
    
    p /= tgamma(k);
    
    return (1 - p);
}

static double igf(double s, double z){
    
    int i;
    double sc=0, sum=1, nom=1, denom=1;
    
    if(z < 0)
        return 0;
    
    sc = 1 / s;
    
    sc *= pow(z, s);
    sc *= exp(-z);
    
    for(i=0;i<200;i++){
        nom *= z;
        s++;
        denom *= s;
        sum += nom / denom;
    }
    
    return sum * sc;
}

int isNumeric (const char *s){
    
    char *p;
    
    if(s == NULL || *s == '\0' || isspace(*s))
        return 0;
    strtod(s, &p);
    return *p == '\0';
}

void lineTerminator(char *line){
    
    int i;
    
    for(i=0;line[i]!=0;i++){
        if(line[i] == '\n' || line[i] == '\r')
            line[i] = '\0';
    }
}

void printInfo(void){
    
    fprintf(stderr,"\n----------------------------------------------------------------------");
    fprintf(stderr,"\nNetwork Enrichment Test\n");
    fprintf(stderr,"Version %s\n", version);
    fprintf(stderr,"\nContact: tuomas.hamala@gmail.com\n");
    fprintf(stderr,"\nRequired parameters:\n");
    fprintf(stderr,"\t-m --modules [file] file defining modules\n");
    fprintf(stderr,"\t-i --in [file] file containing candidate genes or file containing all genes with their test scores\n");
    fprintf(stderr,"\t-t --test [0,1] the choise of test: [0] candidate genes [1] all genes\n\n");
    fprintf(stderr,"Optional paraments:\n");
    fprintf(stderr,"\t-n --neutral [file] file containing neutral test scores\n");
    fprintf(stderr,"\t-c --cycles [int] number of permutation cycles\n");
    fprintf(stderr,"\t-s --side [0,1] are permutation p-values counted against higher [0] over lower [1] estimates (default 0)\n\n");
    fprintf(stderr,"Usage examples:\n");
    fprintf(stderr,"\t./NET -m example.modules.txt -i example.candidates.txt -t 0 > candidate.output.txt\n");
    fprintf(stderr,"\t./NET -m example.modules.txt -i example.allgenes.txt -t 1 -c 10000 > all.output.txt\n");
    fprintf(stderr,"\t./NET -m example.modules.txt -i example.allgenes.txt -n example.neutral.txt -t 1 -c 1e6 > all.output.txt\n");
    fprintf(stderr,"----------------------------------------------------------------------\n\n");
}
