#ifndef RANDOMWALKWITHRESTARTWITHJOINT
#define RANDOMWALKWITHRESTARTWITHJOINT

#include <string>
#include <fstream>
#include <time.h>

#include "walks/walk.hpp" 
#include "api/datatype.hpp"

/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program.
 */
 
class RandomWalkwithRestartwithJoint : public RandomWalk {

public:
    wid_t R;
    hid_t L;

public:

    void initializeRW(wid_t _R, hid_t _L){
        R = _R;
        L = _L;
    }

    hid_t updateByWalk(WalkDataType walk, wid_t walkid, bid_t exec_block, eid_t *&beg_pos, vid_t *&csr, WalkManager &walk_manager ,std::vector<bool> &used_csr, std::vector<bool> &used_csr_v,std::unordered_map<unsigned int, std::vector<int> > &cache){
            //get current time in microsecond as seed to compute rand_r
            tid_t threadid = omp_get_thread_num();
            WalkDataType nowwalk = walk;
            vid_t sourId = walk_manager.getSourceId(nowwalk);
            vid_t curId = walk_manager.getCurrentId(nowwalk) + blocks[exec_block];
            vid_t dstId = curId;
            hid_t hop = walk_manager.getHop(nowwalk);
            // unsigned seed = (unsigned)std::chrono::high_resolution_clock::now().time_since_epoch().count();
            unsigned seed = walk+curId+hop+(unsigned)time(NULL);
            while (( (dstId >= blocks[exec_block] && dstId < blocks[exec_block+1]) || cache.find(dstId)!=cache.end() ) && hop < L ){
            //while (dstId >= blocks[exec_block] && dstId < blocks[exec_block+1] ){
                updateInfo(sourId, dstId, threadid, hop);
                if(cache.find(dstId)!=cache.end() && !(dstId >= blocks[exec_block] && dstId < blocks[exec_block+1])){
                    bool use_cache = false;
                    for(size_t i=0;i<cache[dstId].size();i++){
                        if(cache[dstId][i] != -1){
                            vid_t tmp=dstId;
                            dstId = cache[dstId][i];
                            use_cache = true;
                            cache[tmp][i] = -1;
                            break;
                        }
                    }
                    if(use_cache){
                        hop++;
                        nowwalk++;
                        continue;
                    }else{
                        //cache.erase(dstId);//多线程会出问题
                    }
                    break;
                }

                vid_t dstIdp = dstId - blocks[exec_block];
                eid_t outd = beg_pos[dstIdp+1] - beg_pos[dstIdp];
                // //IO利用率
                // for(vid_t i=0; i<outd; i++){
                //     used_csr[beg_pos[dstIdp]-beg_pos[0]+i] = true;
                // }
                // used_csr_v[dstIdp] = true;
                if ((dstId >= blocks[exec_block] && dstId < blocks[exec_block+1]) && outd > 0 && (float)rand_r(&seed)/RAND_MAX > 0 ){
                //if (outd > 0 && (float)rand_r(&seed)/RAND_MAX > 0 ){
                    eid_t pos = beg_pos[dstIdp] - beg_pos[0] + ((eid_t)rand_r(&seed))%outd;
                    dstId = csr[pos];
                }else{
                    dstId = sourId;
                }
                hop++;
                nowwalk++;
                if(hop%L == L-1) break;
            }
            if( hop%L != L-1 ){
                bid_t p = getblock( dstId );
                if(p>=nblocks) return hop;
                walk_manager.moveWalk(nowwalk, p, threadid, dstId - blocks[p]);
                walk_manager.setMinStep( p, hop );
                walk_manager.ismodified[p] = true;
            }
            return hop;
    }
};

#endif