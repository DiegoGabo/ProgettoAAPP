#include <bitset> 
#include <iostream>

using namespace std;

/*A Nucleotite: A, C, G or T*/
class Nucleotide
{
	private:
		std::bitset<2> nucleotide;
	
	public:
	 	Nucleotide(char);
		string toString();
		bool equal(Nucleotide);
		int getBit(int);
};
