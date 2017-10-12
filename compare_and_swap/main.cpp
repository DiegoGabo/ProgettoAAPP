#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include "structures.hpp"
#include "HashTable.hpp"
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#define L 10

using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
	char ch;
	int k_lenght = 4;
	std::string file_name;

	//clock initialized for calculating execution time 
    struct timeval start, end;
   	gettimeofday(&start, NULL);

	/*program options: the user can insert the lenght of k or ask help*/
	po::options_description desc;
    
    desc.add_options()
        ("h, help", "Shows description of the options")
		("f, file", po::value<std::string>(&file_name)->default_value("../dna_sequences/DNA_prova.txt"), "Set the name of the file in which there is the dna sequence; default value ../dna_sequences/DNA_prova.txt.")
        ("k, k_lenght", po::value<int>(&k_lenght)->default_value(4), "Set the lenght of k; default value 4.");

	po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
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

	/*cout << "Sequence read:\n";
	for(Nucleotide nucleotide : dna_sequence)
		cout << nucleotide.toString();*/

	HashTable hashTable(k_lenght);	//initialize the hash table

	/*cycle in which the k-mer are added to the hash table and so counted*/
	#pragma omp parallel for simd
	for(int i=0; i<dna_sequence.size()-k_lenght; i++)
	{
		std::vector<Nucleotide> k_mer;
		for(int j=i; j<i+k_lenght; j++)
		{
			k_mer.push_back(dna_sequence[j]);
		}
		hashTable.incrementValue(k_mer);
	}
	//salvo su disco
	{
		std::string result_name = "result";
		result_name.append(std::to_string(hashTable.getNum()));
		std::ofstream outfile(result_name);
		hashTable.order(0, pow(2, L));
		cout << "\nscrivo su " << result_name;
		outfile << hashTable.toString();
		outfile.close();
		hashTable.incrementNum();
	}
	
	std::vector<fstream> results;
	std::vector<int> key(hashTable.getNum());
	std::vector<long int> count(hashTable.getNum());
	//unisci file
	for(int i=0; i<hashTable.getNum(); i++)
	{
		std::string result_name = "result";
		result_name.append(std::to_string(i));
		fstream result_stream(result_name, fstream::in);
		//results.pushBack(result_stream);
	}

	bool finish = false;
	
	while(!finish)
	{
		for(int i=0; i<hashTable.getNum(); i++)
		{
			if(!(results[i] >> key[i] >> count[i]))
			{
				key[i]=-1;
			}
		}
		finish = true;
		for(int i=0; i<hashTable.getNum(); i++)
		{
			if(key[i] != -1)
				finish = false;
		}
	}
	for(int i=0; i<hashTable.getNum(); i++)
		cout << "\nKey " << std::to_string(key[i]) << "\tcount " << std::to_string(count[i]);
	

	/*now the extecution time is calcolated and then printed*/
	gettimeofday(&end, NULL);
    float executionTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;	
    cout << "\nExecution time : "<< executionTime << "\n";
}

