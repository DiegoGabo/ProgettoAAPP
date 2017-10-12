#include <Eigen/Dense>
#include <vector>  
#include <string>

using namespace Eigen;

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
		HashTable(int); 
		void incrementValue(std::vector<Nucleotide>);
		std::string toString();
		void order(int, int);
		int getNum();
		void incrementNum();
		void flush();
};

