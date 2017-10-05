#include "structures.hpp"

using namespace std;

Nucleotide::Nucleotide()
{
}


/*constructor that given a char (A, C, G or T) create the Nucleotide*/
Nucleotide::Nucleotide(char value)
{
	switch(value)
	{
		case ('A'):
		    break;
		case ('a'):
		    break;
		 
		case ('C'):
			nucleotide.set(0);
		    break;
		case ('c'):
			nucleotide.set(0);
		    break;
	 
		case ('G'):
			nucleotide.set(1);
			break;
		case ('g'):
			nucleotide.set(1);
			break;

		case ('T'):
			nucleotide.set();
		    break;
		case ('t'):
			nucleotide.set();
		    break;
	}
}

/*check if two nucleotides are equal*/
bool Nucleotide::equal(Nucleotide n2)
{
	if(nucleotide[0]==n2.getBit(0) && nucleotide[1]==n2.getBit(1))
		return true;
	return false;
}

/*obtain the bit of the nucleotide in the position passed as a parameter */
int Nucleotide::getBit(int pos)
{
	return nucleotide[pos];
}

/*return the string that describes the nucleotide*/
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

