#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include "structures.hpp"
#include "HashTable.hpp"
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#define MAX_FILE 100

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
        ("h, help", "Shows description of the options")
		("l, lenght", po::value<int>(&L)->default_value(10), "Set the lenght of the hash table; default value 10.")
		("t, thread", po::value<int>(&thread)->default_value(24), "Set the number of threads")
		("f, file", po::value<std::string>(&file_name)->default_value("../dna_sequences/DNA_prova.txt"), "Set the name of the file in which there is the dna sequence; default value ../dna_sequences/DNA_prova.txt.")
        ("k, k_lenght", po::value<int>(&k_lenght)->default_value(4), "Set the lenght of k; default value 4.");

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
	#pragma omp parallel for simd
	for(int i=0; i<=dna_sequence.size()-k_lenght; i++)
	{
		std::vector<Nucleotide> k_mer;
		for(int j=i; j<i+k_lenght; j++)
		{
			k_mer.push_back(dna_sequence[j]);
		}
		hashTable.incrementValue(k_mer);
	}
	
	//the remaining entries of the table are saved on disk
	std::string result_name = "result";
	result_name.append(std::to_string(hashTable.getNum()));
	std::ofstream outfile(result_name);
	hashTable.order(0, pow(2, L));
	outfile << hashTable.toString();
	outfile.close();
	hashTable.incrementNum();
	
	
	ifstream results[MAX_FILE];	//array in which there are the file written

	std::vector<int> key(hashTable.getNum());
	std::vector<long int> count(hashTable.getNum());

	std::ofstream out_final_file("final_result"); //file in which the result will be written
	
	//all the file previously written are now opened
	for(int i=0; i<hashTable.getNum(); i++)
	{
		std::string result_name = "result";
		result_name.append(std::to_string(i));
		results[i].open(result_name);
	}
	
	//the first line of all file is read
	for(int i=0; i<hashTable.getNum(); i++)
		if(!(results[i] >> key[i] >> count[i]))
			key[i]=-1;		
	
	bool finish = false;

	//while alll file are completely read the cicle sum the values of key that are equal
	while(!finish)
	{
		int min = INT_MAX;
		//the minimum key is found
		for(int i=0; i<hashTable.getNum(); i++)
			if(key[i] < min && key[i] != -1)
				min = key[i];	
		
		int count_min = 0;
		//the count of the minimun key in icremented taking into account all files
		for(int i=0; i<hashTable.getNum(); i++)
			if(key[i] == min)
				do
				{
					count_min += count[i];
					if(!(results[i] >> key[i] >> count[i]))
						key[i]=-1;
				}
				while(key[i] != -1 && key[i]==min);

		finish = true;
		//cycle that check if all files are completely read
		for(int i=0; i<hashTable.getNum(); i++)
			if(key[i] != -1)
				finish = false;

		//translate the key
		std::string mer;
		min--;
		out_final_file << "\n" << min;
		for(int i=0;i<k_lenght;i++){
			int bit_zero = min%2;
			min = min >> 1;
			int bit_one = min%2;
			min = min >> 1;
			char ch = Nucleotide::getChar(bit_zero,bit_one);
			mer.push_back(ch);
		}
		
		out_final_file << "\t" << mer << "\t" << std::to_string(count_min);
	}
	
	/*now the extecution time is calcolated and then printed*/
	gettimeofday(&end, NULL);
    float executionTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;	
    cout << "\nExecution time : "<< executionTime << "\n";
}

