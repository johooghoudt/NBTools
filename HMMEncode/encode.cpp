#include "nbhmm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <clx/table.h>
#include <clx/tokenizer.h>
#include <libgen.h>

#define DIRMODE (S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)

struct Config
{
  const std::string dir;
  int seed;
  double df; // IW's variance parameter
  double dS; // means variance parameter
  double kappa;
  double alpha;
  double gamma;
  int state;
  int iteration;
  
  Config(const char *dir) : dir(dir) {}
};

void ReadConfig(Config &config)
{
  std::cin >> config.seed;
  std::cin >> config.df;
  std::cin >> config.dS;
  std::cin >> config.kappa;
  std::cin >> config.alpha;
  std::cin >> config.gamma;
  std::cin >> config.state;
  std::cin >> config.iteration;
}

dgematrix ReadDataFile(const char *filename, bool hasHeader = false)
{
  std::ifstream in(filename);
  std::string buf;

  
  if(hasHeader){
    getline(in, buf);
  }


  clx::char_separator<char> sep(',');
  clx::table<std::string> table(in, clx::create_tokenizer<std::string>(sep));
  dgematrix matrix((int) table.size(), table[0].size());
  for (size_t i = 0; i < table.size(); i++){
    for(unsigned int j = 0; j < table[i].size(); j++){
      std::stringstream convertor(table[i][j]);
      convertor >> matrix(i, j);
    }
  }

  return matrix;
}

void WriteResults(const NBShdpHmm &hmm, const std::vector< std::vector<int> > &outputs, const std::vector<std::string> &filenames)
{
  char filename[BUFSIZ];

  for(unsigned int i = 0; i < hmm.G.size(); i++){
    sprintf(filename, "mu-%d.dco", i);
    hmm.G[i].Mu.write(filename);
    sprintf(filename, "sig-%d.dge", i);
    hmm.G[i].Sig.write(filename);
  }

  for(unsigned int i = 0; i < outputs.size(); i++){
    std::vector<int> out = outputs[i];
    std::ofstream outfile(filenames[i].c_str());
    for(std::vector<int>::const_iterator it = out.begin(); 
	it != out.end(); ++it){
      outfile << *it << std::endl;
    }
  }
}

void MakeAndChangeDirectory(const char *dir)
{
  if(mkdir(dir, DIRMODE) != 0){
    std::cerr << "can't make directory " << dir << std::endl;
    exit(1);
  }
  if(chdir(dir) != 0){
    std::cerr << "can't change directory " << dir << std::endl;
    exit(1);
  }
}

void usage(const char *program)
{
  std::cout << "usage: " << program << " DIRECTORY INPUT..." << std::endl;
}

int main(int argc, char *argv[])
{
  if(argc < 2){
    usage(argv[0]);
    return 1;
  }
  
  Config config(argv[1]);
  ReadConfig(config);
  
  srand(config.seed);
  srandom(config.seed);

  std::vector<dgematrix> observations;
  std::vector<std::string> filenames;
  int input_count(argc - 2), min_dim(100);
  
  for(int i = 2; i < argc; i++){
    dgematrix data = ReadDataFile(argv[i]);
    observations.push_back(data);
    min_dim = std::min(min_dim, (int)data.n);
    std::string filename = basename(argv[i]);
    filenames.push_back(filename);
  }
  

  NBShdpHmm estimate;
  estimate.resize(config.state, min_dim);
  for(int i = 0; i < config.state; i++){
    estimate.G[i].hp_A.identity();
    estimate.G[i].hp_A = config.df * estimate.G[i].hp_A;
    estimate.G[i].hp_S.identity();
    estimate.G[i].hp_S = config.dS * estimate.G[i].hp_S;
    estimate.hp_kappa = config.kappa;
    estimate.hp_alpha = config.alpha;
    estimate.hp_gamma = config.gamma;
    estimate.G[i].Mu = MultiGaussSampler(estimate.G[i].hp_m,estimate.G[i].hp_S);
  }

  std::vector< std::vector<int> > outputs;
  outputs.resize(input_count);

  MakeAndChangeDirectory(config.dir.c_str());
  for(int i = 0; i < config.iteration; i++){
    for(int j = 0; j < input_count; j++){
      dgematrix backward;
      backward = BackwardFiltering(estimate, observations[j]);
      outputs[j] = ForwardSampling(estimate, backward, observations[j]);
    }
    estimate.Update_shdp_multi(observations, outputs);

    char dir[BUFSIZ];
    sprintf(dir, "iteration-%d", i);
    MakeAndChangeDirectory(dir);
    WriteResults(estimate, outputs, filenames);
    chdir("..");
  }

  return 0;
}
