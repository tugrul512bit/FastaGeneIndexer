#include "FastaGeneIndexer.h"
#include "lib/CpuBenchmarker.h"
#include <tuple>
#include <algorithm>


struct Sequence
{
	size_t size;
	size_t id;
};

// compares 4 FASTA files (of chromosomes each 3.2GB), uses 5.5GB Video-RAM and ~1.5GB RAM, initializes in 180seconds, computes in 70 seconds, on fx8150 cpu
// finds all common sequences (that are found on all files given by "files" vector)
// example: if file 1 has ACGTACGT and if all other files have same, then result contains ACGTACGT
//          if file 1 has ACCCCTT and if same for all other files, then result also contains ACCCCTT
//          complexity:
//             linear time for different length sequences (simply skips inequal length comparisons)
//             O(N*M) when there are N number of k-length sequences in array-1 and M number of k-length sequences in array-2
//             O(k) for string - string comparison (assuming == operator uses SSE/AVX)
int test()
{
    try
    {
        std::vector<std::string> files = {	"./data/homo_sapiens_chromosome.fa", // 3.2GB each
        									"./data/homo_sapiens_chromosome.fa",
        									"./data/homo_sapiens_chromosome.fa",
        									"./data/homo_sapiens_chromosome.fa"};
    	const int n = files.size();


    	// { <file1,sequence1> -> common }
    	std::map<std::tuple<int,int>,bool> result;
    	std::map<std::tuple<int,int>,bool> result2;
    	std::vector<Sequence> sortedSeq[n];


    	// total bytes read during test
    	size_t tot = 0;
        bool debug = true;
        FastaGeneIndexer cache[n];

        {
        	CpuBenchmarker bench(0,"Init");
			#pragma omp parallel for
			for(int i=0;i<n;i++)
			{
				cache[i]=FastaGeneIndexer(files[i], debug,true,true);
			}
        }


        {
        	CpuBenchmarker bench(0,"------------------- Finding common sequence(total) -----------------------");

        	// sort
        	for(int i=0;i<n;i++)
        	{
        		sortedSeq[i] = std::vector<Sequence>();
        		for(size_t j=0;j<cache[i].n();j++)
        		{
        			Sequence seq;
        			seq.size=cache[i].getSequenceLength(j);
        			seq.id=j;
        			sortedSeq[i].push_back(seq);
        		}
        		std::sort(sortedSeq[i].begin(),sortedSeq[i].end(),[](const auto & e1, const auto & e2){ return e1.size < e2.size;});
        	}


        	auto intersection = [&](int cacheId1, int cacheId2)
			{
        		std::mutex m;

        		CpuBenchmarker bench(0,"------- Finding common sequence(file vs file) ------------");
				size_t ctr1 = 0;
				size_t ctr2 = 0;

				// intersection
				while( (ctr1 < cache[cacheId1].n()) && (ctr2 < cache[cacheId2].n()))
				{
					// no need to fetch/decode/compare different sized strings

					// if a size = size block found
					if(sortedSeq[cacheId1][ctr1].size == sortedSeq[cacheId2][ctr2].size)
					{
						size_t ctrtmp1x = 0;
						size_t ctrtmp2x = 0;

						// find n1
						while(sortedSeq[cacheId1][ctr1+ctrtmp1x].size == sortedSeq[cacheId2][ctr2].size)
						{
							ctrtmp1x++;
						}
						size_t ctrtmp1block = ctrtmp1x;

						// find n2
						while(sortedSeq[cacheId1][ctr1].size == sortedSeq[cacheId2][ctr2+ctrtmp2x].size)
						{
							ctrtmp2x++;
						}
						size_t ctrtmp2block = ctrtmp2x;


						for(size_t ctrtmp1 = 0; ctrtmp1<ctrtmp1block; ctrtmp1++)
						{
							auto data1 = cache[cacheId1].getSequenceParallel(sortedSeq[cacheId1][ctr1+ctrtmp1].id);
							tot += sortedSeq[cacheId1][ctr1+ctrtmp1].size;
							for(size_t ctrtmp2 = 0; ctrtmp2<ctrtmp2block; ctrtmp2++)
							{

								tot += sortedSeq[cacheId2][ctr2+ctrtmp2].size;
								// since sizes equal, direct comparison doable
								if( data1 == cache[cacheId2].getSequenceParallel(sortedSeq[cacheId2][ctr2+ctrtmp2].id))
								{

									// cache-m's ctr1 seq = cache-n's ctr2 seq
									result[std::tuple<int,int>(cacheId1,ctr1+ctrtmp1)]=true;
								}
							}
						}


						ctr2 += ctrtmp2block;
						ctr1 += ctrtmp1block;

						continue;
					}



					if(sortedSeq[cacheId1][ctr1].size<sortedSeq[cacheId2][ctr2].size)
					{
						ctr1++;
					}
					else
					{
						ctr2++;
					}
				}
			};


           	auto intersectionWithMap = [&](std::map<std::tuple<int,int>,bool> & mapIn,std::map<std::tuple<int,int>,bool> & mapOut, int cacheId)
    			{
           			CpuBenchmarker bench(0,"------- Finding common sequence(map vs file) ------------");

           			mapOut.clear();

    				size_t ctr1 = 0;
    				size_t ctr2 = 0;

    				std::vector<Sequence> previousVec;
    				for(const auto & e:mapIn)
    				{
    					Sequence seq;
    					seq.id =std::get<1>(e.first);
    					seq.size=cache[std::get<0>(e.first)].getSequenceLength(std::get<1>(e.first));
    					previousVec.push_back(seq);
    				}
    				std::sort(previousVec.begin(), previousVec.end(),[](const Sequence & e1,const Sequence & e2){return e1.size<e2.size;});

    				// intersection
    				while( (ctr1 < previousVec.size()) && (ctr2 < cache[cacheId].n()))
    				{
    					// no need to fetch/decode/compare different sized strings

    					// if a size = size block found
    					if(previousVec[ctr1].size == sortedSeq[cacheId][ctr2].size)
    					{
    						size_t ctrtmp1x = 0;
    						size_t ctrtmp2x = 0;

    						// find n1
    						while(previousVec[ctr1+ctrtmp1x].size == sortedSeq[cacheId][ctr2].size)
    						{
    							ctrtmp1x++;
    						}
    						size_t ctrtmp1block = ctrtmp1x;

    						// find n2
    						while(previousVec[ctr1].size == sortedSeq[cacheId][ctr2+ctrtmp2x].size)
    						{
    							ctrtmp2x++;
    						}
    						size_t ctrtmp2block = ctrtmp2x;

    						for(size_t ctrtmp1 = 0; ctrtmp1<ctrtmp1block; ctrtmp1++)
    						{
    							auto data1 = cache[0].getSequenceParallel(previousVec[ctr1+ctrtmp1].id);
    							tot += previousVec[ctr1+ctrtmp1].size;
    							for(size_t ctrtmp2 = 0; ctrtmp2<ctrtmp2block; ctrtmp2++)
    							{
    								tot += sortedSeq[cacheId][ctr2+ctrtmp2].size;
    								// since sizes equal, direct comparison doable
    								if(data1 == cache[cacheId].getSequenceParallel(sortedSeq[cacheId][ctr2+ctrtmp2].id))
    								{
    									// cache-m's ctr1 seq = cache-n's ctr2 seq
    									mapOut[std::tuple<int,int>(0,ctr1+ctrtmp1)]=true;
    								}
    							}
    						}


    						ctr2 += ctrtmp2block;
    						ctr1 += ctrtmp1block;

    						continue;
    					}



    					if(previousVec[ctr1].size<sortedSeq[cacheId][ctr2].size)
    					{
    						ctr1++;
    					}
    					else
    					{
    						ctr2++;
    					}
    				}
    			};


        	intersection(0,1);
        	intersectionWithMap(result,result2,2);
        	intersectionWithMap(result2,result,3);
        }



        std::cout<<"num sequence found in intersection: "<<result.size()<<std::endl;
        for(const auto & e:result)
        {
        	std::cout<<"sequence:"<<std::endl;
        	if(cache[std::get<0>(e.first)].getSequenceLength(std::get<1>(e.first))<100)
        		std::cout<<(cache[std::get<0>(e.first)].getSequence(std::get<1>(e.first)))<<std::endl;
        	else
        		std::cout<<(cache[std::get<0>(e.first)].getSequence(std::get<1>(e.first),0,100)+std::string("..."))<<std::endl;
        	std::cout<<"---------"<<std::endl;
        }

        std::cout<<"total bytes read after init: "<<tot<<std::endl;
    }
    catch(std::exception & e)
    {
        std::cout<< e.what() <<std::endl;
    }
    return 0;
}
