#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <vector>

using std::endl;
using std::cin;
using std::cout;

void write_csv(std::string filename, std::vector<int32_t> values) {
    std::ofstream outfile(filename);
    
    for (int i = 0; i < values.size(); ++i) {
    if (i < values.size() - 1) {
        outfile << values[i] << "\n";
        }
    else
        outfile << values[i];
    }

    outfile.close();
}

int main() {
    std::default_random_engine generator;
    std::uniform_int_distribution<int32_t> distribution(1, INT32_MAX - 1);
    int n_samples;
    std::string filename;
    cout << "Seeds to generate: " << endl;
    cin >> n_samples;
    cout << "Filepath and name: " << endl;
    cin >> filename;

    std::vector<int32_t> seeds;

    for (int i = 0; i < n_samples; ++i) {
        int32_t gen = distribution(generator);
        seeds.emplace_back(gen);
    }

    write_csv(filename, seeds);

    return 0;
}