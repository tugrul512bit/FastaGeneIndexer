#include"examples/FastaGeneIndexer.h"
#include"CpuBenchmarker.h"
#include<omp.h>
#include<random>
#include<atomic>
int main(int argC, char** argV)
{
	std::string fileName = std::string(argV[argC-1]);

	try
	{
		FastaGeneIndexer fasta(fileName,true);
		std::cout<<fasta.getDescriptor(0)<<std::endl;
		std::cout<<fasta.getSequence(0)<<std::endl;

		for(int j=0;j<100;j++)
		{
			CpuBenchmarker bench;
			const int n = fasta.n();
			std::atomic<size_t> total;
			total.store(0);


			#pragma omp parallel for
			for(int i=0;i<n;i++)
			{
				//for(int k=0;k<5;k++)
				{
					std::string tmp = fasta.getSequence(i);
					total.fetch_add(tmp.length());
				}
			}
			std::cout<<total.load()<<" bytes read"<<std::endl;
		}
	}
	catch (std::exception& ex)
	{
		std::cout << ex.what();

	}
	return 0;
}
