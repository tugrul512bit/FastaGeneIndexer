#include"examples/FastaGeneIndexer.h"

int main(int argC, char** argV)
{
	std::string fileName = std::string(argV[argC-1]);
	try
	{
		bool debug = true;
		FastaGeneIndexer fasta(fileName,debug);
		std::cout<<fasta.getDescriptor(0)<<std::endl;
		std::cout<<fasta.getSequence(0)<<std::endl;		
	}
	catch (std::exception& ex)
	{
		std::cout << ex.what();
	}
	return 0;
}
