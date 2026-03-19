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
 
class RandomWalkwithProb : public RandomWalk {

public:  

    hid_t updateByWalk(WalkDataType walk, wid_t walkid, bid_t exec_block, eid_t *&beg_pos, vid_t *&csr, WalkManager &walk_manager ,
        std::vector<bool> &used_csr, std::vector<bool> &used_csr_v,std::unordered_map<unsigned int, std::vector<int> > &cache, eid_t *&beg_static, vid_t *&csr_static){ //, VertexDataType* vertex_value){
        tid_t threadid = omp_get_thread_num();
        WalkDataType nowWalk = walk;
        vid_t sourId = walk_manager.getSourceId(nowWalk);
        vid_t dstId = walk_manager.getCurrentId(nowWalk) + blocks[exec_block];
        hid_t hop = walk_manager.getHop(nowWalk);
        unsigned seed = (unsigned)(walkid+dstId+hop+(unsigned)time(NULL));
        while (( (dstId >= blocks[exec_block] && dstId < blocks[exec_block+1]) || cache.find(dstId)!=cache.end() ) && hop < L ){
        //while (dstId >= blocks[exec_block] && dstId < blocks[exec_block+1] ){
            updateInfo(sourId, dstId, threadid, hop);
            
            auto it = cache.find(dstId);
            if (it != cache.end() && (dstId < blocks[exec_block] || dstId >= blocks[exec_block + 1])) {
                auto& cached_walks = it->second;
                // 【多线程安全改造】：原子性地获取并自增索引
                // 返回的是自增前的值，相当于原子的 size_t tmp = cached_walks[0]++;
                size_t tmp = __sync_fetch_and_add(&cached_walks[0], 1);
                
                if (tmp >= cached_walks.size()) {
                    break; 
                }

                dstId = cached_walks[tmp];
                hop++;
                nowWalk++;
                continue;
            }

            vid_t dstIdp = dstId - blocks[exec_block];
            eid_t outd = beg_pos[dstIdp+1] - beg_pos[dstIdp];
            // //IO利用率
            // for(vid_t i=0; i<outd; i++){
            //     used_csr[beg_pos[dstIdp]-beg_pos[0]+i] = true;
            // }
            // used_csr_v[dstIdp] = true;
            //if (outd > 0 && (float)rand_r(&seed)/RAND_MAX > 0 ){
            if ((dstId >= blocks[exec_block] && dstId < blocks[exec_block+1]) && outd > 0 && (float)rand_r(&seed)/RAND_MAX > 0 ){
                eid_t pos = beg_pos[dstIdp] - beg_pos[0] + ((eid_t)rand_r(&seed))%outd;
                dstId = csr[pos];
            }else{
                return hop;
            }
            hop++;
            nowWalk++;
        }
        // if( hop < L ){
            bid_t p = getblock( dstId );
            if(p>=nblocks) return hop;
            if(p == exec_block){
                logstream(LOG_DEBUG) << "dstId = " << dstId << ", p " << p << std::endl;
                assert(false);
            }
            walk_manager.moveWalk(nowWalk, p, threadid, dstId - blocks[p]);
            walk_manager.setMinStep( p, hop );
            walk_manager.ismodified[p] = true;
        // }
        return hop;
    }

};

#endif