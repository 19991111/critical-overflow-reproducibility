#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <csignal>
#include <sys/resource.h>
#include <unistd.h>
using u64=std::uint64_t;using ld=long double;
static const int QS[]={64,128,256,512},RS[]={128,128,256,1024};
u64 mix(u64 z){z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);}
u64 word(u64 seed,u64 i){return mix(seed+0x9e3779b97f4a7c15ULL*(i+1));}
struct WaveLaw{std::array<ld,50>x,y;std::array<ld,49>a,b;int shift;};
struct Wave{
 std::array<std::array<ld,100>,257> tab{};
 explicit Wave(const std::string&path){std::ifstream f(path);int a,b,c,d;f>>a>>b>>c>>d;if(a!=257||b!=50||c!=-24||d!=24)throw std::runtime_error("wave shape");std::string s;for(auto&r:tab)for(auto&v:r){f>>s;v=strtold(s.c_str(),nullptr);}if(!f)throw std::runtime_error("wave bytes");}
 WaveLaw law(u64 mass,int q,int cap=0)const{WaveLaw z;u64 p=(mass%q)*256;int j=int(p/q);ld t=ld(p%q)/q;z.shift=int(mass/q)-cap;for(int k=0;k<50;k++){z.x[k]=(1-t)*tab[j][k]+t*tab[j+1][k];z.y[k]=(1-t)*tab[j][k+50]+t*tab[j+1][k+50];}for(int k=0;k<49;k++){z.a[k]=z.x[k]-z.x[k+1];z.b[k]=z.y[k]-z.y[k+1];}return z;}
 static ld H(const WaveLaw&u,const WaveLaw&v){ld z=0;for(int k=0;k<49;k++){int r=k+u.shift-v.shift;ld tail=r<0?1:r>=50?0:v.y[r];z+=u.a[k]*tail;}return z;}
 static ld sig(const WaveLaw&u){ld z=0;for(int k=0;k<49;k++)z+=(u.a[k]+u.b[k])*std::max(k-24+u.shift,0);return z/2;}
};
struct State{int q,m,c;std::vector<int>a;u64 M[2]={0,0},Z=0;explicit State(int q_):q(q_),m(q_/2),c(q_/2),a(2*q_,0){}
 void edge(int l,int r){int win=a[l]<=a[r]?l:r;if(a[win]>=c)Z++;a[win]++;M[win/q]++;}
 void check(u64 n)const{u64 total=0,spill=0,mass[2]={0,0};for(int i=0;i<2*q;i++){total+=a[i];mass[i/q]+=a[i];spill+=std::max(a[i]-c,0);}if(total!=n||spill!=Z||mass[0]!=M[0]||mass[1]!=M[1])throw std::runtime_error("mass/spill invariant");}
};
struct Stats{ld square=0,profile_square=0,pair_square=0,maxspill=0;int feedbacks=0,windows=0;long long Dcrit=0;u64 Zcrit=0;};
std::string hist(const std::vector<int>&a){std::map<int,int> h;for(int x:a)h[x]++;std::ostringstream o;bool first=true;for(auto[x,n]:h){if(!first)o<<';';first=false;o<<x<<':'<<n;}return o.str();}
std::string pairs(const State&s,int d){std::map<std::pair<int,int>,int> h;for(int i=0;i<s.m;i++)h[{s.a[d*s.q+i],s.a[(1-d)*s.q+s.m+i]}]++;std::ostringstream o;bool first=true;for(auto&[p,n]:h){if(!first)o<<';';first=false;o<<p.first<<':'<<p.second<<':'<<n;}return o.str();}
long long Hnum(const std::vector<int>&L,const std::vector<int>&R){long long z=0;std::size_t i=0;for(int r:R){while(i<L.size()&&L[i]<=r)i++;z+=i;}return z;}
ld W1(const std::vector<int>&values,const std::array<ld,50>&tails,int shift){std::vector<int> keys=values;for(int k=0;k<50;k++)keys.push_back(k-24+shift);std::sort(keys.begin(),keys.end());keys.erase(std::unique(keys.begin(),keys.end()),keys.end());ld z=0;std::size_t i=0;for(std::size_t j=1;j<keys.size();j++){int left=keys[j-1];while(i<values.size()&&values[i]<=left)i++;int k=left+1-shift+24;ld tail=k<0?1:k>=50?0:tails[k];ld Femp=ld(i)/values.size();z+=std::abs(Femp-(1-tail))*(keys[j]-left);}return z;}
void observe(const State&s,const Wave&w,Stats&t,std::ofstream&o,int rid,int lambda,int kernel,u64 n,int fid,int hid){
 std::array<std::vector<int>,4> a;for(int v=0;v<4;v++){a[v].assign(s.a.begin()+v*s.m,s.a.begin()+(v+1)*s.m);std::sort(a[v].begin(),a[v].end());}
 long long prod=Hnum(a[0],a[3])-Hnum(a[2],a[1]),fixed=0;for(int i=0;i<s.m;i++){fixed+=(s.a[i]<=s.a[s.q+s.m+i]);fixed-=(s.a[s.q+i]<=s.a[s.m+i]);}
 long long kn=kernel==0?fixed:prod,kd=kernel==0?s.m:s.m*s.m;ld gp=ld(prod)/(s.m*s.m),gk=ld(kn)/kd;auto u=w.law(s.M[0],s.q),v=w.law(s.M[1],s.q);ld gw=Wave::H(u,v)-Wave::H(v,u);ld prof=W1(a[0],u.x,u.shift)+W1(a[1],u.y,u.shift)+W1(a[2],v.x,v.shift)+W1(a[3],v.y,v.shift);
 auto cu=w.law(s.M[0],s.q,s.c),cv=w.law(s.M[1],s.q,s.c);ld sp=Wave::sig(cu)+Wave::sig(cv),err=ld(s.Z)/s.q-sp;
 if(fid>=0){t.square+=(gk-gw)*(gk-gw);t.profile_square+=(gp-gw)*(gp-gw);t.pair_square+=(gk-gp)*(gk-gp);t.feedbacks++;}
 if(hid>=0){t.maxspill=std::max(t.maxspill,std::abs(err));t.windows++;}
 o<<rid<<','<<s.q<<','<<lambda<<','<<kernel<<','<<n<<','<<fid<<','<<hid<<','<<s.M[0]<<','<<s.M[1]<<','<<s.Z<<','<<kn<<','<<kd<<','<<prod<<','<<s.m*s.m<<','<<gw<<','<<sp<<','<<prof;
 for(auto&x:a){o<<','<<hist(x);}o<<','<<(kernel==0?pairs(s,0):"")<<','<<(kernel==0?pairs(s,1):"")<<'\n';
}
void put32(std::ofstream&f,std::uint32_t x){for(int k=0;k<4;k++)f.put(char((x>>(8*k))&255));}
void put64(std::ofstream&f,u64 x){for(int k=0;k<8;k++)f.put(char((x>>(8*k))&255));}
void config(std::ofstream&f,const State&s,int rid,int lambda,int kernel,u64 n){put32(f,s.q);put32(f,lambda);put32(f,rid);put32(f,kernel);put64(f,n);for(int v:s.a)put32(f,v);}
void selftest(const Wave&w){
 if(std::numeric_limits<ld>::digits<64||std::numeric_limits<ld>::radix!=2||std::numeric_limits<ld>::round_style!=std::round_to_nearest)throw std::runtime_error("long double contract");
 State s(4);std::vector<int> capped(8,0);u64 spill=0;for(int n=0;n<128;n++){int l=(n/4)%2*4+n%2,r=(n/8)%2*4+2+(n/2)%2;int win=capped[l]<=capped[r]?l:r;if(capped[win]>=s.c)spill++;else capped[win]++;s.edge(l,r);s.check(n+1);for(int i=0;i<8;i++)if(capped[i]!=std::min(s.a[i],s.c))throw std::runtime_error("cap coupling");if(spill!=s.Z)throw std::runtime_error("cap spill");}
 for(int lambda:{0,1,2})for(int kernel:{0,1}){std::array<int,8> inc{};int internal[2]={0,0},cross[2]={0,0};for(int z=0;z<4;z++)for(int d=0;d<2;d++)for(int i=0;i<2;i++)for(int j=0;j<2;j++){bool cr=z<lambda;int lc=d,rc=cr?1-d:d,jj=cr&&kernel==0?i:j;inc[lc*4+i]++;inc[rc*4+2+jj]++;if(cr)cross[d]++;else internal[d]++;}for(int x:inc)if(x!=8)throw std::runtime_error("uniform endpoint");for(int d=0;d<2;d++)if(internal[d]!=(4-lambda)*4||cross[d]!=lambda*4)throw std::runtime_error("type probability");}
 std::vector<int>L={0,1,1,4},R={1,1,2,3};long long naive=0;for(int l:L)for(int r:R)naive+=l<=r;if(Hnum(L,R)!=naive)throw std::runtime_error("integer comparison");auto a=w.law(0,64);if(Wave::H(a,a)-Wave::H(a,a)!=0)throw std::runtime_error("wave symmetry");std::cout<<"{\"status\":\"selftest_passed\",\"random_trajectories\":0,\"long_double_digits\":"<<std::numeric_limits<ld>::digits<<"}\n";
}
int main(int argc,char**argv){try{
 if(argc<3)throw std::runtime_error("usage: executable --selftest wave OR --run wave seeds output");Wave wave(argv[2]);if(std::string(argv[1])=="--selftest"){selftest(wave);return 0;}if(argc!=5||std::string(argv[1])!="--run")throw std::runtime_error("arguments");
 selftest(wave);struct rlimit lim{360,361};setrlimit(RLIMIT_CPU,&lim);struct rlimit mem{1024ULL*1024*1024,1024ULL*1024*1024};setrlimit(RLIMIT_AS,&mem);alarm(420);
 std::ifstream sf(argv[3]);std::vector<u64>seeds;int id;std::string hex;while(sf>>id>>hex){if(id!=int(seeds.size()))throw std::runtime_error("seed order");seeds.push_back(std::stoull(hex,nullptr,16));}if(seeds.size()!=1024)throw std::runtime_error("seed count");
 std::filesystem::path out(argv[4]);if(std::filesystem::exists(out/"outcomes.csv"))throw std::runtime_error("refuse overwrite");std::ofstream result(out/"outcomes.csv"),obs(out/"observations_part00.csv"),cfg(out/"critical_configs.bin",std::ios::binary),progress(out/"progress.csv");if(!result||!obs||!cfg||!progress)throw std::runtime_error("output open");int obs_part=0;u64 obs_closed_bytes=0;
 result<<std::setprecision(21);obs<<std::setprecision(21);result<<"replicate,q,lambda,kernel,g_rmse,profile_g_rmse,pair_g_rmse,max_spill_error,Dcrit_integer,Zcrit,n_completed,feedback_count,window_count,driver_hash\n";obs<<"replicate,q,lambda,kernel,n,feedback_id,window_id,M0,M1,Z,g_kernel_num,g_kernel_den,g_product_num,g_product_den,g_wave,S_wave,profile_W1,L0,R0,L1,R1,pairs01,pairs10\n";progress<<"q,lambda,replicate,arrivals,driver_hash\n";cfg.write("FQCFG1\0",7);
 for(int cell=0;cell<4;cell++){int q=QS[cell],m=q/2,bits=0;while((1<<bits)<q)bits++;int mbits=bits-1;u64 nstar=u64(q)*q,end=nstar+2*q;std::map<u64,std::pair<int,int>> checks;
  for(int j=0;j<16;j++){u64 n=nstar/2+u64((q/4-1)*j/15)*2*q+u64(2*q*j/16);checks[n].first=j+1;}
  for(int j=0;j<5;j++){u64 n=static_cast<long long>(nstar)+static_cast<long long>(j-2)*q;checks[n].second=j+1;}
  for(int lambda:{0,1,2})for(int rid=0;rid<RS[cell];rid++){
   State fixed(q),product(q);Stats ft,pt;u64 digest=1469598103934665603ULL;auto next=checks.begin();
   for(u64 n=1;n<=end;n++){u64 z=word(seeds[rid],2*(n-1)),v=word(seeds[rid],2*(n-1)+1);digest=(digest^z)*1099511628211ULL;digest=(digest^v)*1099511628211ULL;bool cross=int(z&(q-1))<lambda;int d=(z>>bits)&1,i=v&(m-1),j=(v>>mbits)&(m-1);int l=d*q+i,r=(cross?1-d:d)*q+m+j,rf=(cross?1-d:d)*q+m+(cross?i:j);fixed.edge(l,rf);product.edge(l,r);
    if(next!=checks.end()&&n==next->first){fixed.check(n);product.check(n);if(lambda==0&&(fixed.a!=product.a||fixed.Z!=product.Z))throw std::runtime_error("lambda0 mismatch");int fid=next->second.first-1,hid=next->second.second-1;observe(fixed,wave,ft,obs,rid,lambda,0,n,fid,hid);observe(product,wave,pt,obs,rid,lambda,1,n,fid,hid);next++;}
    if(n==nstar){ft.Dcrit=static_cast<long long>(fixed.M[0])-static_cast<long long>(fixed.M[1]);pt.Dcrit=static_cast<long long>(product.M[0])-static_cast<long long>(product.M[1]);ft.Zcrit=fixed.Z;pt.Zcrit=product.Z;config(cfg,fixed,rid,lambda,0,n);config(cfg,product,rid,lambda,1,n);}
   }
   for(int k=0;k<2;k++){auto&t=k?pt:ft;if(t.feedbacks!=16||t.windows!=5)throw std::runtime_error("checkpoint count");result<<rid<<','<<q<<','<<lambda<<','<<k<<','<<sqrtl(t.square/16)<<','<<sqrtl(t.profile_square/16)<<','<<sqrtl(t.pair_square/16)<<','<<t.maxspill<<','<<t.Dcrit<<','<<t.Zcrit<<','<<end<<','<<t.feedbacks<<','<<t.windows<<','<<digest<<'\n';}
   result.flush();obs.flush();cfg.flush();progress<<q<<','<<lambda<<','<<rid<<','<<end<<','<<digest<<'\n';progress.flush();if(u64(result.tellp())+obs_closed_bytes+u64(obs.tellp())+u64(cfg.tellp())>256ULL*1024*1024)throw std::runtime_error("output budget");
   if(u64(obs.tellp())>32ULL*1024*1024){obs_closed_bytes+=u64(obs.tellp());obs.close();std::ostringstream name;name<<"observations_part"<<std::setw(2)<<std::setfill('0')<<++obs_part<<".csv";obs.open(out/name.str());obs<<std::setprecision(21)<<"replicate,q,lambda,kernel,n,feedback_id,window_id,M0,M1,Z,g_kernel_num,g_kernel_den,g_product_num,g_product_den,g_wave,S_wave,profile_W1,L0,R0,L1,R1,pairs01,pairs10\n";if(!obs)throw std::runtime_error("observation chunk open");}
  }
  std::cout<<"{\"q_completed\":"<<q<<",\"replicates_per_lambda\":"<<RS[cell]<<"}"<<std::endl;
 }
 std::cout<<"{\"status\":\"planned_batch_completed\",\"paired_trajectories\":4608,\"SDE_reference_samples\":0}"<<std::endl;return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
