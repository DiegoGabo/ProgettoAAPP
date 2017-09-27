#include <iostream>
#include <fstream>
#include <bitset> 

using namespace std;


class Nucleotide
{
	private:
		 std::bitset<2> nucleotide;
	
	public:
	 	Nucleotide(char);
		char toString();
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

char Nucleotide::toString()
{
	if(nucleotide[0]==0)
	{
		if(nucleotide[1]==0)
			return 'A';
		return 'C';
	}

	if(nucleotide[1]==0)
		return 'G';
	return 'T';
}

int main()
{
	char ch;
	fstream dna_sequence("DNA_prova.txt", fstream::in);
	while (dna_sequence >> noskipws >> ch) 
	{
		Nucleotide nucleotide(ch);
		cout << nucleotide.toString();
	}
}
