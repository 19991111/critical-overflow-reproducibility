#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using ld=long double;using u64=std::uint64_t;
constexpr int N=16384,R=32768,NG=4097,NS=16385;constexpr ld H=1.0L/N,SQH=1.0L/128,PI=0xc.90fdaa22168c235p-2L,LN2=0xb.17217f7d1cf79acp-4L;
std::array<ld,24> odd;std::array<ld,35> invfact;
void init(){for(int j=0;j<24;j++)odd[j]=1.0L/(2*j+1);ld fact=1;invfact[0]=1;for(int j=1;j<35;j++){fact*=j;invfact[j]=1/fact;}}
ld atanseries(ld t){ld t2=t*t,p=t,s=t;for(int j=1;j<24;j++){p*=t2;s+=p*odd[j];}return 2*s;}
ld ownlog(ld x){if(!(x>0))throw std::runtime_error("log domain");if(x>=.5L&&x<=2)return atanseries((x-1)/(x+1));int e;ld m=std::frexp(x,&e)*2;--e;return atanseries((m-1)/(m+1))+e*LN2;}
ld ownsqrt(ld x){if(x==0)return 0;if(x<0)throw std::runtime_error("sqrt domain");int e;ld m=std::frexp(x,&e);if(e%2){m*=2;--e;}ld s=1;for(int j=0;j<7;j++)s=(s+m/s)/2;return std::ldexp(s,e/2);}
std::pair<ld,ld> normals(u64 a,u64 b){constexpr u64 mask=(1ULL<<52)-1;ld u=(ld(a&mask)+.5L)*0x1p-52L,v=(ld(b&mask)+.5L)*0x1p-52L;ld radius=ownsqrt(-2*ownlog(u));int q=int(4*v);ld angle=(4*v-q)*(PI/2),a2=angle*angle,cs=invfact[32],sn=invfact[33];for(int j=15;j>=0;j--){cs=invfact[2*j]-a2*cs;sn=invfact[2*j+1]-a2*sn;}sn*=angle;ld c,s;switch(q){case 0:c=cs;s=sn;break;case 1:c=-sn;s=cs;break;case 2:c=-cs;s=-sn;break;default:c=sn;s=-cs;}return{radius*c,radius*s};}
u64 mix(u64 z){z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);}
u64 word(u64 seed,u64 j){return mix(seed+0x9e3779b97f4a7c15ULL*(j+1));}
struct Tables{
 std::array<ld,NG>g;std::array<ld,NS>s;
 explicit Tables(const char*path){std::ifstream f(path);int a,b,c;f>>a>>b>>c;if(a!=NG||b!=NS||c!=512)throw std::runtime_error("table shape");std::string x;for(auto&v:g){f>>x;v=strtold(x.c_str(),nullptr);}for(auto&v:s){f>>x;v=strtold(x.c_str(),nullptr);}if(!f)throw std::runtime_error("table read");if(g[0]!=0||s.back()!=16||s[0]>=.632L)throw std::runtime_error("table endpoints");for(int i=1;i<NG;i++)if(g[i]>g[i-1]||512*(g[i]-g[i-1])< -2)throw std::runtime_error("g slopes");for(int i=1;i<NS;i++)if(s[i]<s[i-1]||512*(s[i]-s[i-1])>.5001L)throw std::runtime_error("S slopes");}
 ld drift(ld x)const{ld a=std::abs(x),z;if(a>=8)z=g.back();else{ld p=a*512;int j=int(p);z=g[j]+(p-j)*(g[j+1]-g[j]);}return x<0?-z:z;}
 ld readout(ld x)const{ld a=std::abs(x);if(a>=32)return a/2;ld p=a*512;int j=int(p);return s[j]+(p-j)*(s[j+1]-s[j]);}
};
void selftest(const Tables&t){if(std::numeric_limits<ld>::digits<64||std::numeric_limits<ld>::radix!=2||std::numeric_limits<ld>::round_style!=std::round_to_nearest)throw std::runtime_error("arithmetic contract");ld worst=0;for(int e=-53;e<=10;e++)for(int j=0;j<257;j++){ld x=std::ldexp(1+ld(j)/257,e);worst=std::max(worst,std::abs(ownlog(x)-std::log(x)));if(std::abs(ownsqrt(x)-std::sqrt(x))>1e-12L)throw std::runtime_error("sqrt test");}if(worst>1e-12L)throw std::runtime_error("log test");constexpr u64 mask=(1ULL<<52)-1;for(u64 a:std::array<u64,5>{0,1,mask/2,mask-1,mask})for(u64 b:std::array<u64,6>{0,1,mask/4,mask/2,mask-1,mask}){auto z=normals(a,b);ld u=(ld(a)+.5L)*0x1p-52L,v=(ld(b)+.5L)*0x1p-52L,r=std::sqrt(-2*std::log(u));ld c=r*std::cos(2*PI*v),s=r*std::sin(2*PI*v);if(std::abs(z.first-c)>1e-11L||std::abs(z.second-s)>1e-11L)throw std::runtime_error("normal test");}if(t.drift(0)!=0||t.readout(32)!=16)throw std::runtime_error("lookup test");std::cout<<"{\"status\":\"selftest_passed\",\"random_paths\":0,\"long_double_digits\":"<<std::numeric_limits<ld>::digits<<",\"log_check_max_difference\":"<<std::setprecision(20)<<worst<<"}\n";}
volatile ld checksum=0;
void calibrate(const Tables&t){std::clock_t begin=std::clock();constexpr u64 mask=(1ULL<<52)-1;ld z=0;for(u64 i=0;i<262144;i++){auto n=normals((i*0x9e3779b97f4a7c15ULL)&mask,(i*0xd1b54a32d192ed03ULL+17)&mask);ld x=ld(int(i%8192)-4096)/512;z+=(x+H*t.drift(x)+SQH*n.first)+(x+2*H*t.drift(x+.125L)+SQH*n.first)+(x+H*t.drift(x-.25L)+SQH*n.second)+(x+2*H*t.drift(x+.375L)+SQH*n.second);}checksum=z;ld cpu=ld(std::clock()-begin)/CLOCKS_PER_SEC;std::cout<<std::setprecision(20)<<"{\"calibration_pairs\":262144,\"CPU_seconds\":"<<cpu<<",\"main_CPU_prediction_with_factor1_5\":"<<cpu*1024*1.5L<<",\"random_paths\":0,\"checksum\":"<<z<<"}\n";}
void sample(const Tables&t,const char*seedpath,const char*outpath){std::ifstream f(seedpath);std::vector<u64> seeds;int id;std::string s;while(f>>id>>s){if(id!=int(seeds.size()))throw std::runtime_error("seed order");seeds.push_back(std::stoull(s,nullptr,16));}if(seeds.size()!=R)throw std::runtime_error("seed count");std::ifstream check(outpath);if(check.good())throw std::runtime_error("refuse overwrite");std::ofstream out(outpath),progress(std::string(outpath)+".progress.csv");if(!out||!progress)throw std::runtime_error("output open");out<<"replicate,lambda,D_hex,S_hex,driver_hash\n";progress<<"replicate,CPU_seconds\n";std::clock_t start=std::clock();for(int rid=0;rid<R;rid++){ld a=0,b=0;u64 digest=1469598103934665603ULL;for(int j=0;j<N/2;j++){u64 u=word(seeds[rid],2*j),v=word(seeds[rid],2*j+1);digest=(digest^u)*1099511628211ULL;digest=(digest^v)*1099511628211ULL;auto z=normals(u,v);a+=H*t.drift(a)+SQH*z.first;b+=2*H*t.drift(b)+SQH*z.first;a+=H*t.drift(a)+SQH*z.second;b+=2*H*t.drift(b)+SQH*z.second;}if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a)>1155||std::abs(b)>1155)throw std::runtime_error("state bound");out<<rid<<",1,"<<std::hexfloat<<a<<','<<t.readout(a)<<','<<std::defaultfloat<<digest<<'\n';out<<rid<<",2,"<<std::hexfloat<<b<<','<<t.readout(b)<<','<<std::defaultfloat<<digest<<'\n';if((rid+1)%256==0){out.flush();progress<<rid+1<<','<<ld(std::clock()-start)/CLOCKS_PER_SEC<<'\n';progress.flush();if((rid+1)%4096==0)std::cout<<"{\"reference_replicates_completed\":"<<rid+1<<"}"<<std::endl;}}
 std::cout<<"{\"status\":\"all_reference_paths_completed\",\"replicates_per_lambda\":"<<R<<",\"steps\":"<<N<<",\"CPU_seconds\":"<<ld(std::clock()-start)/CLOCKS_PER_SEC<<"}"<<std::endl;
}
ld capital(const std::vector<ld>&x,int k,bool upper){ld m=ld(k)/100000,s=0,c=0;for(ld v:x){ld factor=1+(upper?m-v:v-m)/8;ld y=ownlog(factor)-c;ld t=s+y;c=(t-s)-y;s=t;}return s;}
std::pair<int,int> interval(const std::vector<ld>&x){ld threshold=ownlog(1000)+.000001L;int l=0,u=300000,L,U;if(capital(x,0,false)<threshold)L=0;else{while(u-l>1){int k=(l+u)/2;if(capital(x,k,false)>=threshold)l=k;else u=k;}L=l;}l=0;u=300000;if(capital(x,u,true)<threshold)U=300000;else{while(u-l>1){int k=(l+u)/2;if(capital(x,k,true)>=threshold)u=k;else l=k;}U=u;}return{L,U};}
void means(const char*path){std::ifstream f(path);std::string line;std::getline(f,line);std::vector<ld>a,b;while(std::getline(f,line)){std::stringstream ss(line);std::vector<std::string>v;std::string s;while(std::getline(ss,s,','))v.push_back(s);if(v.size()!=5)throw std::runtime_error("reference row");ld x=strtold(v[3].c_str(),nullptr);(v[1]=="1"?a:b).push_back(std::min(x,3.0L));}if(a.size()!=R||b.size()!=R)throw std::runtime_error("reference completeness");auto x=interval(a),y=interval(b);std::cout<<"{\"lambda1\":{\"lower_grid\":"<<x.first<<",\"upper_grid\":"<<x.second<<"},\"lambda2\":{\"lower_grid\":"<<y.first<<",\"upper_grid\":"<<y.second<<"},\"grid_denominator\":100000,\"log_capital_error_allowance\":0.000001,\"threshold\":1000}"<<std::endl;}
int main(int argc,char**argv){try{init();if(argc==3&&std::string(argv[1])=="--means"){means(argv[2]);return 0;}if(argc<3)throw std::runtime_error("arguments");Tables t(argv[2]);std::string mode=argv[1];if(mode=="--selftest")selftest(t);else if(mode=="--calibrate")calibrate(t);else if(mode=="--sample"&&argc==5)sample(t,argv[3],argv[4]);else throw std::runtime_error("mode");return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
