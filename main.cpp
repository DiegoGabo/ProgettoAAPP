#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include "structures.hpp"
#include "HashTable.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
	char ch;
	int k_lenght = 4;

	/*program options: the user can insert the lenght of k or ask help*/
	po::options_description desc;
    
    desc.add_options()
        ("help, h", "Shows description of the options")
        ("k-lenght, k", po::value<int>(&k_lenght)->default_value(4), "Set the lenght of k; default value 4.");

	po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if(vm.count("help"))
    {
        cout << desc;
        return 0;
    }
    
	/*the file in which there is the dna sequence is opened*/
	fstream dna_sequence_file("DNA_prova.txt", fstream::in);
	std::vector<Nucleotide> dna_sequence;	

	while (dna_sequence_file >> noskipws >> ch) 
	{
		Nucleotide nucleotide(ch);
		dna_sequence.push_back(nucleotide);
	}

	for(Nucleotide nucleotide : dna_sequence)
	{
		cout << nucleotide.toString();
	}

	HashTable hashTable(k_lenght);

	for(int i=0; i<dna_sequence.size()-k_lenght; i++)
	{
		std::vector<Nucleotide> k_mer;
		for(int j=i; j<i+k_lenght; j++)
		{
			k_mer.push_back(dna_sequence[j]);
		}
		hashTable.incrementValue(k_mer);
	}

	cout << "\nConteggio:\n" << hashTable.toString();
}


