#include <iostream>
#include <fstream>
#include <bitset>
#include <vector>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

class Nucleotide
{
	private:
		 std::bitset<2> nucleotide;

	public:
	 	Nucleotide(char);
		string toString();
		int getBit(int);
		bool equal(Nucleotide);
};

Nucleotide::Nucleotide(char value)
{
	switch(value)
	{
		case ('A'):
		    break;

		case ('C'):
			nucleotide.set(0);
		    break;

		case ('G'):
			nucleotide.set(1);
			break;

		case ('T'):
			nucleotide.set();
		    break;
	}
}

bool Nucleotide::equal(Nucleotide n2)
{
	if(nucleotide[0]==n2.getBit(0) && nucleotide[1]==n2.getBit(1))
		return true;
	return false;
}

int Nucleotide::getBit(int pos)
{
	return nucleotide[pos];
}

string Nucleotide::toString()
{
	if(nucleotide[0]==0)
	{
		if(nucleotide[1]==0)
			return "A";
		return "G";
	}

	if(nucleotide[1]==0)
		return "C";
	return "T";
}


int main(int argc, char *argv[])
{
	char ch;
	int k_lenght = 4;

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

	cout << "\n\nK_mers:\n\n";

	for(int i=0; i<dna_sequence.size()-k_lenght; i++)
	{
		std::vector<Nucleotide> k_mer;
		for(int j=i; j<i+k_lenght; j++)
		{
			k_mer.push_back(dna_sequence[j]);
			cout << dna_sequence[j].toString();
		}
		cout << "\n";
	}
}


