#include "conf.h"

int verifyFileExistance (const char *fname);
tspcfg loadConfigurationFile(char *fname);
void parseCfgLine(tspcfg *config, char *line);
int blankLine(char *line);
void printUsage();

tspcfg argParse(int argc, char **argv) {

    //verify if an help message was requested
    for(int i=1; i<argc; i++){
        if( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 ){
            printUsage();
            exit(0);
        }
    }
    
    //verify if another file location is provided via command line argument
    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--config-file") == 0) {
            if (argc > i+1) {
                char *fname = argv[i+1];
                if (verifyFileExistance(fname)) {
                    printf("Configuration loaded from: %s\n", fname);
                    return loadConfigurationFile(fname);
                } else {
                    printf("Configuration file %s not found.\n", fname);
                    exit(1);
                }    
            }else{
                printf("No filename provided after option -c\n");
                exit(1);
            }
        }
    }

    //set up default file location list
    char filelist[3][50] = {
        "gims.conf"
        "$HOME/.conf/gims.conf",
        "/etc/gims/gims.conf",
    };

    for (int i = 0; i < 3; i++) {
        if (verifyFileExistance(filelist[i])) {
            printf("Loading configuration from: %s\n", filelist[i]);
            return loadConfigurationFile(filelist[i]);
        }
    }

    printf("Could not locate any configuration file to load.");
    exit(1);
}

void printUsage() {
    printf(
        "tsp-evo - Travelling Salesman Problem evolutionary solver\n"
        "\n"
        "Usage: tsp-evo [-h|--help] [-c|--config-file <configuration file>]\n"
        "\n"
        "Options:\n"
        "  -h, --help                    display this help message and exit\n"
        "  -c FILE, --config-file=FILE   load configurations from FILE\n"
    );
}

int verifyFileExistance (const char *fname) {
  struct stat buffer;
  return (stat (fname, &buffer) == 0); 
}

tspcfg loadConfigurationFile(char *fname) {
    
    tspcfg config;
    FILE *f = fopen(fname, "r");

    char *line = (char *)malloc(256 * sizeof(char));
    size_t size = sizeof(char) * 256;

    int nbytes;

    while ((nbytes = getline(&line, &size, f)) != -1) {
        if(line[0] != '#' && !blankLine(line)){
            if(line[nbytes-1] == '\n')
                line[nbytes-1] = '\0';
            parseCfgLine(&config, line);
        }
    }

    free(line);
    fclose(f);
    return config;
}

void parseCfgLine(tspcfg *config, char *line) {

    int i=0;
    while(line[i] != ' ')
        i++;

    char *value = line+i+1;
    line[i] = '\0';

    if (strcmp(line, "nthreads") == 0) {
        config->nthreads = atoi(value);
    }
    else if (strcmp(line, "popsize") == 0) {
        config->popsize = atoi(value);
    }
    else if (strcmp(line, "npops") == 0) {
        config->npops = atoi(value);
    }
    else if (strcmp(line, "ngens") == 0) {
        config->ngens = atoi(value);
    }
    else if (strcmp(line, "maxgrad0count") == 0) {
        config->maxgrad0count = atoi(value);
    }
    else if (strcmp(line, "swap_mutation_rate") == 0) {
        config->swap_mutation_rate = atof(value);
    }
    else if (strcmp(line, "inversion_mutation_rate") == 0) {
        config->inversion_mutation_rate = atof(value);
    }
    else if (strcmp(line, "individual_replacement_rate") == 0) {
        config->individual_replacement_rate = atof(value);
    }
    else if (strcmp(line, "elitesize") == 0) {
        config->elitesize = atoi(value);
    }
    else if (strcmp(line, "tournament_size") == 0) {
        config->tournament_size = atoi(value);
    }
    else if (strcmp(line, "crossover_rate") == 0) {
        config->crossover_rate = atof(value);
    }
    else if (strcmp(line, "migrations") == 0) {
        config->migrations = atoi(value);
    }
    else if (strcmp(line, "population_migration_rate") == 0) {
        config->population_migration_rate = atof(value);
    }
    else if (strcmp(line, "individual_migration_rate") == 0) {
        config->individual_migration_rate = atof(value);
    }
    else{
        printf("unknown keyword in configuration file: %s\n", line);
        exit(1);
    }
}

int blankLine(char *line){
    int i=0;
    while( line[i] == '\t' || line[i] == ' ' )
        i++;
    if( line[i] != '\n' )
        return 0;
    return 1;
}
