#include"FastaGeneIndexer.h"

int main(int argC, char** argV)
{
	try
	{
		bool debug = true;
		FastaGeneIndexer fasta("./data/influenza.fna",debug);
		std::cout<<fasta.getDescriptor(0)<<std::endl;
		std::cout<<fasta.getSequence(0)<<std::endl;	
		
		fasta.initDescriptorIndexMapping();
		std::string sequence = fasta.getSequenceByDescriptor("gi|58576|gb|X52226|Influenza A virus (A/FPV/Rostock/34(H7N1)) gene for neuraminidase, genomic RNA");		
		std::cout<< sequence <<std::endl;
	}
	catch (std::exception& ex)
	{
		std::cout << ex.what();
	}
	return 0;
}
