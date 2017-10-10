#include <Eigen/Dense>
#include <vector>  
#include <string>
#include <atomic>
#include <cstdlib>

#define K 4
#define L 8

using namespace Eigen;

/*a single entry of the hash table*/
class HashEntry
{
	private:
		Nucleotide key[K];
		int count;
	public:
		HashEntry() noexcept;
		bool isEmpty();
		int getC();
		Nucleotide* getK();
		void setK(Nucleotide key[]);
		void setC(int);
		std::string toString();
		int* getPC();
};

class HashTable
{
	private:
		std::atomic<HashEntry> table[L];
		MatrixXf matrix, inverse;
		int f(Nucleotide []);
		int reprobe(int);
	public:
		HashTable(); 
		void incrementValue(Nucleotide n[]);
		std::string toString();
};
