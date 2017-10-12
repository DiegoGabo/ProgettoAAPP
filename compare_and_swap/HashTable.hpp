#include <Eigen/Dense>
#include <vector>  
#include <string>

using namespace Eigen;

/*class used to implement the hash table in which are counted the k-mers*/
class HashTable
{
	private:
		int k;
		int m;
		int num;
		bool stop;
		std::vector<int> table_hash; 
		std::vector<int> table_count; //count (32bit) ..ottimizzazione 2?
		MatrixXf matrix, inverse;
		int f(int);
		int key_vector(std::vector<Nucleotide>);
		int reprobe(int);
		int partition_hash(int, int);
		void swap(int, int);	
		
	public:
		HashTable(int, int); 
		void incrementValue(std::vector<Nucleotide>);
		std::string toString();
		int getNum();
		void incrementNum();
		void flush();
		void order(int, int);
};

