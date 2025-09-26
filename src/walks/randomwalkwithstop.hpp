#ifndef RANDOMWALKWITHSTOP
#define RANDOMWALKWITHSTOP

#include <string>
#include <fstream>
#include <time.h>

#include "walks/walk.hpp" 
#include "api/datatype.hpp"

/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program.
 */
 
class RandomWalkwithStop : public RandomWalk {

public:  

    hid_t updateByWalk(WalkDataType walk, wid_t walkid, bid_t exec_block, eid_t *&beg_pos, vid_t *&csr, WalkManager &walk_manager ,std::unordered_map<unsigned int, std::vector<int> > &cache){ //, VertexDataType* vertex_value){
        tid_t threadid = omp_get_thread_num();
        WalkDataType nowWalk = walk;
        vid_t sourId = walk_manager.getSourceId(nowWalk);
        vid_t dstId = walk_manager.getCurrentId(nowWalk) + blocks[exec_block];
        hid_t hop = walk_manager.getHop(nowWalk);
        unsigned seed = (unsigned)(walkid+dstId+hop+(unsigned)time(NULL));
        while (( (dstId >= blocks[exec_block] && dstId < blocks[exec_block+1]) || cache.find(dstId)!=cache.end() ) && hop < L ){
        //while (( (dstId >= blocks[exec_block] && dstId < blocks[exec_block+1]) ) && hop < L ){
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
                    nowWalk++;
                    continue;
                }
                break;
            }
            vid_t dstIdp = dstId - blocks[exec_block];
            eid_t outd = beg_pos[dstIdp+1] - beg_pos[dstIdp];
            if ((dstId >= blocks[exec_block] && dstId < blocks[exec_block+1]) && outd > 0 && (float)rand_r(&seed)/RAND_MAX > 0 ){
                eid_t pos = beg_pos[dstIdp] - beg_pos[0] + ((eid_t)rand_r(&seed))%outd;
                dstId = csr[pos];
            }else{
                if(outd == 0) {
                    return hop;
                }
                break;
            }
            hop++;
            nowWalk++;
        }
        if( hop < L ){
            bid_t p = getblock( dstId );
            if(p>=nblocks) return hop;
            walk_manager.moveWalk(nowWalk, p, threadid, dstId - blocks[p]);
            walk_manager.setMinStep( p, hop );
            walk_manager.ismodified[p] = true;
        }
        return hop+1;
    }

};

#endif