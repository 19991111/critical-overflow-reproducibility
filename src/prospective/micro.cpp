#include <algorithm>
#include <array>
#include <cfenv>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <sched.h>
using u64=std::uint64_t;
constexpr int Q=512,M=256,C=256,R=4096;
u64 mix(u64 z){z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);}
u64 word(u64 seed,u64 i){return mix(seed+0x9e3779b97f4a7c15ULL*(i+1));}
struct State {
 int q,m,c;std::vector<int>a;u64 mass[2]={0,0},spill=0;
 explicit State(int q_):q(q_),m(q_/2),c(q_/2),a(2*q_,0){}
 void edge(int l,int r){int w=a[l]<=a[r]?l:r;if(a[w]>=c)++spill;++a[w];++mass[w/q];}
 void check(u64 n)const{u64 total=0,z=0,mm[2]={0,0};for(int j=0;j<2*q;++j){total+=a[j];z+=std::max(a[j]-c,0);mm[j/q]+=a[j];}if(total!=n||z!=spill||mm[0]!=mass[0]||mm[1]!=mass[1])throw std::runtime_error("mass/spill invariant");}
};
void environment(){
 cpu_set_t cpus;CPU_ZERO(&cpus);if(sched_getaffinity(0,sizeof(cpus),&cpus))throw std::runtime_error("affinity");int count=CPU_COUNT(&cpus);
 unsigned short cw=0;unsigned int mxcsr=0;
#if defined(__x86_64__) || defined(__i386__)
 asm volatile("fnstcw %0":"=m"(cw));asm volatile("stmxcsr %0":"=m"(mxcsr));
#endif
 if(count!=1||std::fegetround()!=FE_TONEAREST||std::numeric_limits<long double>::radix!=2||std::numeric_limits<long double>::digits<64)throw std::runtime_error("environment contract");
 std::cout<<"{\"status\":\"environment_before_words\",\"affinity_count\":"<<count<<",\"FE_TONEAREST\":true,\"long_double_digits\":"<<std::numeric_limits<long double>::digits<<",\"x87_control_word\":"<<cw<<",\"MXCSR\":"<<mxcsr<<"}"<<std::endl;
}
void selftest(){
 State s(4);std::vector<int> capped(8,0);u64 zs=0;
 for(int n=0;n<128;++n){int l=(n/4)%2*4+n%2,r=(n/8)%2*4+2+(n/2)%2,w=capped[l]<=capped[r]?l:r;if(capped[w]>=s.c)++zs;else ++capped[w];s.edge(l,r);s.check(n+1);for(int j=0;j<8;++j)if(capped[j]!=std::min(s.a[j],s.c))throw std::runtime_error("cap coupling");if(zs!=s.spill)throw std::runtime_error("cap spill");}
 State tie(4);tie.edge(0,6);if(tie.a[0]!=1||tie.a[6]!=0)throw std::runtime_error("global left tie");
 for(int lam:{0,1,2})for(int kernel:{0,1}){std::array<int,8>inc{};int cross[2]={0,0},internal[2]={0,0};for(int z=0;z<4;++z)for(int d=0;d<2;++d)for(int i=0;i<2;++i)for(int j=0;j<2;++j){bool cr=z<lam;int rc=cr?1-d:d,jj=cr&&kernel==0?i:j;++inc[d*4+i];++inc[rc*4+2+jj];if(cr)++cross[d];else ++internal[d];}for(int x:inc)if(x!=8)throw std::runtime_error("endpoint law");for(int d=0;d<2;++d)if(cross[d]!=lam*4||internal[d]!=(4-lam)*4)throw std::runtime_error("edge law");}
 std::cout<<"{\"status\":\"selftest_passed\",\"stochastic_trajectories\":0}"<<std::endl;
}
void put32(std::ofstream&f,std::uint32_t x){for(int k=0;k<4;++k)f.put(char((x>>(8*k))&255));}
void put64(std::ofstream&f,u64 x){for(int k=0;k<8;++k)f.put(char((x>>(8*k))&255));}
void config(std::ofstream&f,const State&s,int rid,int lam,int kernel){put32(f,Q);put32(f,lam);put32(f,rid);put32(f,kernel);put64(f,u64(Q)*Q);for(int x:s.a)put32(f,x);}
int main(int argc,char**argv){try{
 environment();if(argc==2&&std::string(argv[1])=="--selftest"){selftest();return 0;}
 if(argc!=4||std::string(argv[1])!="--run")throw std::runtime_error("usage --selftest OR --run seeds empty-output-directory");
 std::ifstream sf(argv[2]);int id;std::string hex;std::vector<u64> seeds;while(sf>>id>>hex){if(id!=int(seeds.size()))throw std::runtime_error("seed order");seeds.push_back(std::stoull(hex,nullptr,16));}if(!sf.eof()||seeds.size()!=R)throw std::runtime_error("seed count/parse");
 std::filesystem::path out(argv[3]);if(!std::filesystem::exists(out)||!std::filesystem::is_empty(out))throw std::runtime_error("refuse nonempty output");
 std::ofstream rows(out/"outcomes.csv"),progress(out/"progress.csv");if(!rows||!progress)throw std::runtime_error("output open");
 rows<<"replicate,q,lambda,kernel,Dcrit_integer,Zcrit,n_completed,driver_hash\n";progress<<"lambda,replicate,arrivals,driver_hash\n";
 for(int lam:{1,2}){std::ofstream cfg;
  for(int rid=0;rid<R;++rid){if(rid%1024==0){cfg.close();std::string name="terminal_configs_lambda"+std::to_string(lam)+"_part"+std::to_string(rid/1024)+".bin";cfg.open(out/name,std::ios::binary);if(!cfg)throw std::runtime_error("config open");cfg.write("PMCFG1\0",7);}
   State fixed(Q),product(Q);u64 digest=1469598103934665603ULL;
   for(u64 n=0;n<u64(Q)*Q;++n){u64 z=word(seeds[rid],2*n),v=word(seeds[rid],2*n+1);digest=(digest^z)*1099511628211ULL;digest=(digest^v)*1099511628211ULL;bool cr=int(z&(Q-1))<lam;int d=(z>>9)&1,i=v&(M-1),j=(v>>8)&(M-1),l=d*Q+i;fixed.edge(l,(cr?1-d:d)*Q+M+(cr?i:j));product.edge(l,(cr?1-d:d)*Q+M+j);}
   for(int k=0;k<2;++k){const State&s=k?product:fixed;s.check(u64(Q)*Q);config(cfg,s,rid,lam,k);long long d=static_cast<long long>(s.mass[0])-static_cast<long long>(s.mass[1]);rows<<rid<<','<<Q<<','<<lam<<','<<k<<','<<d<<','<<s.spill<<','<<u64(Q)*Q<<','<<digest<<'\n';}
   rows.flush();cfg.flush();progress<<lam<<','<<rid<<','<<u64(Q)*Q<<','<<digest<<'\n';progress.flush();if(!rows||!cfg||!progress)throw std::runtime_error("write failure");
   if((rid+1)%1024==0)std::cout<<"{\"lambda\":"<<lam<<",\"paired_replicates_completed\":"<<rid+1<<"}"<<std::endl;
  }
 }
 std::cout<<"{\"status\":\"complete_once\",\"paired_runs\":8192,\"endpoint_rows\":16384,\"new_SDE_paths\":0}"<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
