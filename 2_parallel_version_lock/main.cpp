#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include "structures.hpp"
#include "HashTable.hpp"
#include <time.h>
#include <sys/time.h>
#include <omp.h>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
	char ch;
	int k_lenght = 4;
	int L = 10;
	int thread = 24;
	std::string file_name;

	//clock initialized for calculating execution time 
    struct timeval start, end;
   	gettimeofday(&start, NULL);

	/*program options: the user can insert the lenght of k or ask help*/
	po::options_description desc;
    
    desc.add_options()
        ("help, h", "Shows description of the options")
		("l, lenght", po::value<int>(&L)->default_value(10), "Set the lenght of the hash table; default value 10.")
		("t, thread", po::value<int>(&thread)->default_value(24), "Set the number of threads")
		("file, f", po::value<std::string>(&file_name)->default_value("../dna_sequences/DNA_prova.txt"), "Set the name of the file in which there is the dna sequence; default value ../dna_sequences/DNA_prova.txt.")
        ("k-lenght, k", po::value<int>(&k_lenght)->default_value(4), "Set the lenght of k; default value 4.");

	po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

	omp_set_num_threads(thread);
    
    if(vm.count("help"))
    {
        cout << desc;
        return 0;
    }
    
	/*the file in which there is the dna sequence is opened*/
	fstream dna_sequence_file(file_name, fstream::in);
	std::vector<Nucleotide> dna_sequence;	

	/*cycle in which the dna_sequence is read and codified as a vector of nucleotides*/
	while (dna_sequence_file >> noskipws >> ch) 
	{
		if(ch=='A' || ch=='a' || ch=='C' || ch=='c' || ch=='G' || ch=='g' || ch=='T' || ch=='t')
		{
			Nucleotide nucleotide(ch);
			dna_sequence.push_back(nucleotide);
		}
	}

	HashTable hashTable(k_lenght, L);	//initialize the hash table

	/*cycle in which the k-mer are added to the hash table and so counted*/
	#pragma omp parallel for
	for(int i=0; i<=dna_sequence.size()-k_lenght; i++)
	{
		std::vector<Nucleotide> k_mer;
		for(int j=i; j<i+k_lenght; j++)
		{
			k_mer.push_back(dna_sequence[j]);
		}
		#pragma omp critical
		hashTable.incrementValue(k_mer);
	}

	//the result is written on the final_result file 
	std::ofstream outfile("final_result");
	outfile << hashTable.toString();

	/*now the extecution time is calcolated and then printed*/
	gettimeofday(&end, NULL);
    float executionTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;	
    cout << "\nExecution time : "<< executionTime << "\n";
}


