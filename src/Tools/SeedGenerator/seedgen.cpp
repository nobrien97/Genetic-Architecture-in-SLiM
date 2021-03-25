#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <getopt.h>

using std::endl;
using std::cout;

#define no_argument 0
#define required_argument 1

void write_csv(std::string filename, std::vector<int32_t> values, std::string header = "Seed") {
    std::ofstream outfile(filename);

    if (header.size()) {
        outfile << header << "\n";
    }
    
    for (int i = 0; i < values.size(); ++i) {
    if (i < values.size() - 1) {
        outfile << values[i] << "\n";
        }
    else
        outfile << values[i];
    }

    outfile.close();
}


void doHelp(char* appname) {
    std::fprintf(stdout,
    "Uniformly Distributed Seed Generator\n"
    "\n"
    "This program generates a .csv of uniformly distributed 32-bit integers for use as RNG seeds.\n"
    "Usage: %s [OPTION]...\n"
    "Example: %s -h\n"
    "\n"
    "-h     Print this help manual.\n"
    "\n"
    "-n N   Generate N random samples.\n"
    "\n"
    "-v     Turn on verbose mode.\n"    
    "-d FILEPATH    Specify a filepath and name for the generated seeds to be saved.\n"
    "Example: -d ~/Desktop/seeds.csv\n"
    "\n",
    appname,
    appname
    );

}

int main(int argc, char* argv[]) {
    

    int optionindex = 0;

    const struct option longopts[] =
    {
        { "destination",    no_argument,        0,  'd' },
        { "nsamples",       required_argument,  0,  'n' },
        { "help",           no_argument,        0,  'h' },
        { "verbose",        no_argument,        0,  'v' },
        {0,0,0,0}
    };

// Initialise variables with defaults if values are not supplied
    std::string filename = "./seeds.csv";
    int n_samples = 10;
    bool debug = false;

    int options; 

    while (options != -1) {

        options = getopt_long(argc, argv, "n:d:hv", longopts, &optionindex);

        switch (options) {
            case 'n':
                n_samples = std::stoi(optarg);
                continue;

            case 'd':
                filename = optarg;
                continue;

            case 'h':
                doHelp(argv[0]);
                return 1;

            case 'v':
                debug = true;
                continue;

            case -1:
                break;
            }
        }


    std::random_device mersseed;
    if(debug) {
        cout << "/dev/random seed for Mersenne Twister: " << mersseed() << "\n"
             << "Number of seeds =  " << n_samples << "\n"
             << "File written to: " << filename << endl;
    }

    std::mt19937 generator(mersseed());
    std::uniform_int_distribution<int32_t> distribution(1, INT32_MAX - 1);

    std::vector<int32_t> seeds;

    for (int i = 0; i < n_samples; ++i) {
        int32_t gen = distribution(generator);
        seeds.emplace_back(gen);
    }

    write_csv(filename, seeds);

    return 0;
}