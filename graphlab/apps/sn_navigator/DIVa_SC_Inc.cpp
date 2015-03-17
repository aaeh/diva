/* This software is distributed under the MIT license (see LICENSE file in the
 * root directory)
 *
 * Copyright (c) 2013 Amira A. Soliman
 */

#include <graphlab.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/ui/metrics_server.hpp>
#include <graphlab/macros_def.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>    // std::reverse
/*
*#################################### Data Structure Part #####################################
*/

size_t snapShotID=0;

///////////////////////////////global functions/////////////////////////////////

int find(std::vector<std::string> vec, std::string key){
    int index =-1;

    for(size_t i=0; i<vec.size(); i++){
       if(vec[i] == key){
          index =i;
          break;
       }
    }
    return index;
}

int count(std::vector<std::string> vec, std::string key){
    int ct =0;

    for(size_t i=0; i<vec.size(); i++){
       if(vec[i] == key){
          ct ++;
       }
    }
    return ct;
}

int count(std::vector<graphlab::vertex_id_type> vec, graphlab::vertex_id_type key){
    int ct =0;

    for(size_t i=0; i<vec.size(); i++){
       if(vec[i] == key){
          ct ++;
       }
    }
    return ct;
}

bool compare(size_t i,size_t j) { return (i>j); }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//gossip message exchanged between social graph nodes
//////////////////////////////////////////////////////////

struct single_message {
    graphlab::vertex_id_type sender_id;
    graphlab::vertex_id_type reciever_id; 
    std::vector<graphlab::vertex_id_type> path;
    std::vector<std::string> att_group;
    std::vector<size_t> freq;
    std::vector<size_t> size;
    

    single_message(): att_group(),freq(),size(), path() { }

   void save(graphlab::oarchive& oarc) const {
                oarc <<sender_id<<reciever_id;
                
                size_t num = att_group.size();
                oarc << num;
                for(std::vector<std::string>::const_iterator iter = att_group.begin();
                       iter != att_group.end(); ++iter){
                               oarc << *iter;
                }

                num = freq.size();
                oarc << num;
                for(std::vector<size_t>::const_iterator iter = freq.begin();
                       iter != freq.end(); ++iter){
                               oarc << *iter;
                }
                num = size.size();
                oarc << num;
                for(std::vector<size_t>::const_iterator iter = size.begin();
                       iter != size.end(); ++iter){
                               oarc << *iter;
                }

                num = path.size();
                oarc << num;
                for(std::vector<graphlab::vertex_id_type>::const_iterator iter = path.begin();
                       iter != path.end(); ++iter){
                               oarc << *iter;
                }

   }
   void load(graphlab::iarchive& iarc) {
          iarc>>sender_id>>reciever_id;
          size_t num,f;
          std::string vk;
          graphlab::vertex_id_type val;

           iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk;
                  att_group.push_back(vk);
               }
        iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> f;
                  freq.push_back(f);
               }
         iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> f;
                  size.push_back(f);
               }

         iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> val;
                  path.push_back(val);
               }
   }
};

// message object exchanged using ivertex program
struct data_message: public graphlab::IS_POD_TYPE {
  std::vector<single_message> msgs;
  double value;

  data_message():  value(1.0), msgs() { }

  //
  // returns a priority to the reception of the message (used by the asynchronous engines)
  //
  double priority() const { return std::fabs(value); }

  //
  // allow message merging across the network
  // GraphLab always merges messages destined to the same vertex.
  //
  data_message& operator+=(const data_message& other) {

          for(size_t i = 0;i < other.msgs.size(); ++i){
               msgs.push_back(other.msgs[i]);
          }

          return *this;
  }
  void save(graphlab::oarchive& oarc) const {

        size_t num_msgs = msgs.size();
        oarc << num_msgs;
        for(size_t k = 0;k < num_msgs; ++k){
                 msgs[k].save(oarc);
        }
  }

  void load(graphlab::iarchive& iarc) {

          size_t num_msgs = 0;
          msgs.clear();

          iarc >> num_msgs;

          for(size_t k = 0;k < num_msgs; ++k){
                  single_message msg;
                  msg.load(iarc);
                  msgs.push_back(msg);
          }
  }
};

//////////////////////////////////////////////////////////
//edge
struct edata {
  graphlab::edge_id_type weight;

  edata():weight(0) {
  }

  explicit edata(graphlab::edge_id_type w): weight(w) {
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << weight;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> weight;
  }

};
////////I//////////////////////////////////////////////////
//vertex
struct vdata {

        // list of direct friends
                
        size_t join_time;
        graphlab::vertex_id_type commID;
        std::vector<graphlab::vertex_id_type> path_to_Diva;
        std::map<graphlab::vertex_id_type, graphlab::vertex_id_type> frnds_commIDs;
        std::map<graphlab::vertex_id_type, std::vector<graphlab::vertex_id_type> > Divas_paths;
        std::vector<graphlab::vertex_id_type> comm_members;
        std::map<graphlab::vertex_id_type, size_t> clusters;
        std::multimap<std::string, std::string>  att_values;
        std::map<graphlab::vertex_id_type, std::vector<std::string> > lCAS;
        std::map<graphlab::vertex_id_type, std::vector<size_t> > lFreq;
        std::map<graphlab::vertex_id_type, std::vector<size_t> > lSize;
        std::map<graphlab::vertex_id_type, std::pair<std::string, double> > commID_CAS;
        

     vdata() : 
               path_to_Diva(), frnds_commIDs(), comm_members(), clusters(),Divas_paths(), att_values(), lCAS(),lFreq(), lSize(), commID_CAS() {}

      explicit vdata(size_t jt, std::multimap<std::string, std::string> att_v) :
          path_to_Diva(), frnds_commIDs(), comm_members(), clusters(), Divas_paths(), att_values(att_v), lCAS(),lFreq(), lSize(), commID_CAS(), join_time(jt) {}


      void find_DIVas(){
           std::vector<size_t> link_count;
           size_t diff =0;
           clusters.clear();
           for(std::map<graphlab::vertex_id_type, graphlab::vertex_id_type>::iterator it= frnds_commIDs.begin();
               it !=frnds_commIDs.end(); ++it){
               if(clusters.find(it->second) == clusters.end())
                  clusters.insert(std::make_pair(it->second,1));
               else
                 clusters[it->second] ++;
           }
          for(std::map<graphlab::vertex_id_type, size_t>::iterator it= clusters.begin();
             it !=clusters.end(); ++it)
              link_count.push_back(it->second);  
                      
          std::sort(link_count.begin(), link_count.end(), compare);  
          for(size_t i=1; i< link_count.size(); i++)
             diff+= link_count[0]-link_count[i];
          if(diff!=0){ //means that there is dominant ID among my friends
             //
            for(std::map<graphlab::vertex_id_type, size_t>::iterator it= clusters.begin();
              it !=clusters.end(); ++it)
              if(it->second == link_count[0]){
                commID = it->first;
                path_to_Diva = Divas_paths[commID];
                break; 
              }
          }
          /*else{ //no dominant ID, select the max ID among communities
             std::map<graphlab::vertex_id_type, size_t>::iterator it= clusters.end();
             it --;
             commID =  it->first;
             path_to_Diva = Divas_paths[commID];
          }*/
      }

      void compute_LCAS (std::map<graphlab::vertex_id_type, std::multimap<std::string, std::string> > new_frnds, std::map<graphlab::vertex_id_type, std::vector<std::string> > new_LCAS, std::map<graphlab::vertex_id_type, std::vector<size_t> > new_lFreq, std::map<graphlab::vertex_id_type, std::vector<size_t> > new_lSize){
          
          std::vector<std::string> att_names;
          att_names.push_back("gndr");
          att_names.push_back("fname");
          att_names.push_back("education");
          att_names.push_back("job");
          att_names.push_back("employer");
          att_names.push_back("interest");
          att_names.push_back("home_cntry");
          att_names.push_back("currloc_cntry");

          for(std::map<graphlab::vertex_id_type, size_t>::iterator clst = clusters.begin(); clst!=clusters.end(); ++clst){
               graphlab::vertex_id_type com = clst->first;
               std::vector<std::string> cas;
               std::vector<size_t> frq, sz;
               std::vector <std::multimap<std::string, std::string> > frnds_prfl;
               for(std::map<graphlab::vertex_id_type, std::multimap<std::string, std::string> >::iterator
                  iter = new_frnds.begin(); iter!= new_frnds.end(); ++iter){
                   if(frnds_commIDs[iter->first] == com)
                     frnds_prfl.push_back(iter->second);
               }

               for(size_t i=0;i<att_names.size()-1; i++){
                 std::string a_n1 = att_names[i];
                 std::vector<std::string> triple;
                 std::vector<std::string> val1, val2, f_val1, f_val2;
                 for(size_t j=i+1; j<att_names.size(); j++){
                    std::string a_n2 = att_names[j];
                    for(size_t c=0; c<frnds_prfl.size(); c++){
                       std::multimap<std::string, std::string> frd_p = frnds_prfl[c];
                       for(std::multimap<std::string, std::string>::iterator it1 = frd_p.equal_range(a_n1).first;
                          it1!= frd_p.equal_range(a_n1).second; ++it1)
                           val1.push_back(it1->second);
                        for(std::multimap<std::string, std::string>::iterator it2 = frd_p.equal_range(a_n2).first;
                          it2!= frd_p.equal_range(a_n2).second; ++it2)
                           val2.push_back(it2->second);

                      size_t tot_tran_size =  val1.size() + val2.size();
                      size_t freq = 1;
                      std::string v1,v2;
                      for(size_t i=0; i<val1.size(); i++){
                         v1 = val1[i];
                         for(size_t j=0; i< val2.size(); j++){
                            v2=val2[j];
                            for(size_t rest_of_frnds = 0;rest_of_frnds < frnds_prfl.size(); rest_of_frnds ++){
                              if(rest_of_frnds != c){
                                 frd_p = frnds_prfl[rest_of_frnds];
                                 for(std::multimap<std::string, std::string>::iterator it1 = frd_p.equal_range(a_n1).first;
                                    it1!= frd_p.equal_range(a_n1).second; ++it1)
                                    f_val1.push_back(it1->second);
                                 for(std::multimap<std::string, std::string>::iterator it2 = frd_p.equal_range(a_n2).first;
                                    it2!= frd_p.equal_range(a_n2).second; ++it2)
                                    f_val2.push_back(it2->second);

                                    if(find(f_val1, v1) >0 && find(f_val2, v2) >0){
                                      size_t c1 = count(f_val1, v1);
                                      size_t c2 = count(f_val1, v2);
                                      if(c1<c2)
                                        freq+=c1;
                                      else freq+=c2;
                                    }
                              }
                            }
                         }
                      }
                     //add a_n1&a_n2 to CAS
                     int in = find(cas, a_n1 + "&" +  a_n2);
                     if(in>0){
                       frq[in] += freq;
                       sz[in] +=tot_tran_size;
                     }
                     else{
                         cas.push_back(a_n1 + "&" +  a_n2);
                         frq.push_back(freq);
                         sz.push_back(tot_tran_size);
                     }
                    }
                 }
             }
               //new_LCASes
               new_LCAS.insert(std::make_pair(com,cas));
               new_lFreq.insert(std::make_pair(com,frq));
               new_lSize.insert(std::make_pair(com,sz));
          }

          
       }

        void save(graphlab::oarchive& oarc) const {
              
                oarc << commID << join_time;

                size_t num = path_to_Diva.size();
                oarc << num;
                for(size_t i = 0;i < num; ++i){
                        oarc << path_to_Diva[i];
                }

                num = frnds_commIDs.size();
                oarc << num;
                for(std::map<graphlab::vertex_id_type, graphlab::vertex_id_type>::const_iterator iter = frnds_commIDs.begin();
                       iter != frnds_commIDs.end(); ++iter){
                        oarc << iter->first << iter->second;
                }
                
                num = comm_members.size();
                oarc << num;
                for(size_t i = 0;i < num; ++i){
                        oarc << comm_members[i];
                }
 
                num = clusters.size();
                oarc << num;
                for(std::map<graphlab::vertex_id_type, size_t>::const_iterator iter = clusters.begin();
                       iter != clusters.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                }

                num = Divas_paths.size();
                oarc << num;
                for(std::map<graphlab::vertex_id_type, std::vector<graphlab::vertex_id_type> >::const_iterator iter = Divas_paths.begin();
                       iter != Divas_paths.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second.size();
                               for(size_t j=0; j< iter->second.size(); j++)
                                   oarc << iter->second[j];
                }

                num = att_values.size();
                oarc << num;
                for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                }
               
                num = lFreq.size();
                oarc << num;
                for(std::map<graphlab::vertex_id_type, std::vector<size_t> >::const_iterator iter=lFreq.begin();
                      iter!= lFreq.end(); iter++){
                      oarc << iter->first;
                      for(size_t i = 0;i < iter->second.size(); ++i){
                         oarc << iter->second[i];
                      } 
                }

                num = lCAS.size();
                oarc << num;
                for(std::map<graphlab::vertex_id_type, std::vector<std::string> >::const_iterator iter=lCAS.begin();
                      iter!= lCAS.end(); iter++){
                      oarc << iter->first;
                      for(size_t i = 0;i < iter->second.size(); ++i){
                         oarc << iter->second[i];
                      } 
                }

                num = lSize.size();
                oarc << num;
                for(std::map<graphlab::vertex_id_type, std::vector<size_t> >::const_iterator iter=lSize.begin();
                      iter!= lSize.end(); iter++){
                      oarc << iter->first;
                      for(size_t i = 0;i < iter->second.size(); ++i){
                         oarc << iter->second[i];
                      }
                }
          //std::multimap<graphlab::vertex_id_type, std::pair<std::string, double> > commID_CAS;
                num = commID_CAS.size();
                oarc << num;
                for(std::map<graphlab::vertex_id_type, std::pair<std::string, double> >::const_iterator iter = commID_CAS.begin();
                       iter != commID_CAS.end(); ++iter){
                               //std::pair<std::string, double> p = iter->second;
                               oarc << iter->first;
                               oarc << iter->second.first << iter->second.second;
                }
             
        }

     void load(graphlab::iarchive& iarc) {
               
               std::string v,vk;
               size_t num = 0, cnt=0;
               graphlab::vertex_id_type nID, cID;
                iarc >> commID >> join_time;

                iarc >> num;
                for(size_t i = 0;i < num; ++i){
                        iarc >> nID; 
                        path_to_Diva.push_back(nID);
                }

                iarc >> num;
                for(size_t i = 0;i < num; ++i){
                        iarc >> nID >> cID;
                        frnds_commIDs.insert(std::make_pair(nID, cID));
                }

                iarc >> num;
                for(size_t i = 0;i < num; ++i){
                        iarc >> nID;
                        comm_members.push_back(nID);
                }

                iarc >> num;
                for(size_t i = 0;i < num; ++i){
                    iarc >> nID >> cnt;
                    clusters.insert(std::make_pair(nID, cnt));
                }
                
                iarc >> num;
                for(size_t i = 0;i < num; ++i){
                   iarc >> nID;
                   std::vector<graphlab::vertex_id_type> p;
                   size_t len;
                   iarc >> len;
                   for(size_t j = 0;j < len; ++j){
                      graphlab::vertex_id_type n;
                      iarc >> n;
                      p.push_back(n);
                   }
                   Divas_paths.insert(std::make_pair(nID, p)); 
                }                

                num = att_values.size();
                iarc >> num;
                for(size_t i = 0;i < num; ++i){
                    iarc >> vk >> v;
                    att_values.insert(std::make_pair(vk, v));
                }

                iarc >> num;
                for(size_t i = 0;i < num; ++i){
                        iarc >> nID;
                        size_t len, v;
                        std::vector<size_t> f;
                        iarc >> len;
                        for(size_t j = 0;j < len; ++j){
                           iarc >> v;
                           f.push_back(v);
                        } 
                        lFreq.insert(std::make_pair(nID, f));
                }

                iarc >> num;
                for(size_t i = 0;i < num; ++i){
                    iarc >> nID;
                    size_t len;
                    std::string v;
                    std::vector<std::string> f;
                    iarc >> len;
                    for(size_t j = 0;j < len; ++j){
                           iarc >> v;
                           f.push_back(v);
                        }
                    lCAS.insert(std::make_pair(nID, f));
                }

                iarc >> num;
                for(size_t i = 0;i < num; ++i){
                    iarc >> nID;
                        size_t len, v;
                        std::vector<size_t> f;
                        iarc >> len;
                        for(size_t j = 0;j < len; ++j){
                           iarc >> v;
                           f.push_back(v);
                        }
                    lSize.insert(std::make_pair(nID, f));
                }

          //std::multimap<graphlab::vertex_id_type, std::pair<std::string, double> > commID_CAS;
                iarc >> num;
                double d;
                for(size_t i = 0;i < num; ++i){
                   iarc >> nID >> vk >> d;
                   commID_CAS.insert(std::make_pair(nID, std::make_pair(vk,d)));
                }
        
        }
 };



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//distributed graph definition
typedef graphlab::distributed_graph<vdata, edata> graph_type;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#################################### Send data Items part #####################################
//

struct data_gatherer_directN {

        std::vector<graphlab::vertex_id_type> frnds_IDs;
        std::map<graphlab::vertex_id_type, std::multimap<std::string, std::string> > frnds_profiles;
        std::map<graphlab::vertex_id_type, graphlab::vertex_id_type> frnds_commIDs;
        std::map<graphlab::vertex_id_type, std::vector<graphlab::vertex_id_type> > path_to_DIVas;
        

       data_gatherer_directN () : frnds_IDs(), frnds_profiles(), frnds_commIDs(), path_to_DIVas() {}
        void save(graphlab::oarchive& oarc) const {
           size_t num = frnds_profiles.size();
                oarc << num;
                
                for(std::map<graphlab::vertex_id_type, std::multimap<std::string, std::string> >::const_iterator iter = frnds_profiles.begin();
                       iter != frnds_profiles.end(); ++iter){
                               oarc << iter->first;
                       
                       for(std::multimap<std::string, std::string>::const_iterator it = iter->second.begin(); it != iter->second.end(); it++){
                               oarc << it->first;
                               oarc << it->second;
                       }
                }
            

            num = frnds_IDs.size();
            for(size_t c=0; c<num; c++)
              oarc <<frnds_IDs[c]; 

            num = frnds_commIDs.size();
            oarc << num;
            for(std::map<graphlab::vertex_id_type, graphlab::vertex_id_type>::const_iterator iter =frnds_commIDs.begin();
               iter != frnds_commIDs.end(); iter++)
                   oarc << iter->first << iter->second;
          
            num = path_to_DIVas.size();
            oarc << num;
            for(std::map<graphlab::vertex_id_type, std::vector<graphlab::vertex_id_type> >::const_iterator iter = path_to_DIVas.begin();
               iter != path_to_DIVas.end(); iter++){
                  oarc << iter->first;
                  oarc << iter->second.size();
                  for(size_t i=0; i< iter->second.size(); i++)
                      oarc << iter->second[i];
            }

        }
        void load(graphlab::iarchive& iarc) {
             size_t num, count;
             graphlab::vertex_id_type f, c;
             frnds_profiles.clear();
             std::string vk, v;
             std::multimap<std::string, std::string> att_values;
             iarc >> count;
             for(int j=0; j<count;j++){
             att_values.clear();
              iarc >> f;
              iarc >> num;
           
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
               }
            frnds_profiles.insert(std::make_pair(f,att_values));
           }

           iarc >> num;
           for(size_t i=0; i<num; i++){
              iarc >> f >> c;
              frnds_commIDs.insert(std::make_pair(f,c));
           }

          iarc >> num;
          for(size_t i=0; i<num; i++){
             iarc >> c;
             size_t plen;
             std::vector<graphlab::vertex_id_type> path;
             iarc >> plen;
             for(size_t j=0; i<plen; j++){
               iarc >>f;
               path.push_back(f);
             }
             path_to_DIVas.insert(std::make_pair(c, path));
          }
        }
   
       data_gatherer_directN& operator+=(const data_gatherer_directN& other) {
              for(std::map<graphlab::vertex_id_type, std::multimap<std::string, std::string> >::const_iterator iter = other.frnds_profiles.begin();
                       iter != other.frnds_profiles.end(); ++iter){
                  frnds_profiles.insert(std::make_pair(iter->first, iter->second));
              }
              for(std::vector<graphlab::vertex_id_type>::const_iterator iter = other.frnds_IDs.begin(); iter != other.frnds_IDs.end(); ++iter){
                 frnds_IDs.push_back(*iter);
              }
              for(std::map<graphlab::vertex_id_type, graphlab::vertex_id_type>::const_iterator itr= other.frnds_commIDs.begin(); itr != other.frnds_commIDs.end(); itr++) 
                  frnds_commIDs.insert(std::make_pair(itr->first, itr->second));

              for(std::map<graphlab::vertex_id_type, std::vector<graphlab::vertex_id_type> >::const_iterator itr= other.path_to_DIVas.begin(); itr !=  other.path_to_DIVas.end(); itr++)
                  path_to_DIVas.insert(std::make_pair(itr->first, itr->second));

              return *this;
         }

};

class DIVa_NodeCentric: public graphlab::ivertex_program<graph_type,data_gatherer_directN,data_message>, public graphlab::IS_POD_TYPE {

public:
 
       std::map<graphlab::vertex_id_type, std::multimap<std::string, std::string> > frnds_profiles;
       std::map<graphlab::vertex_id_type, graphlab::vertex_id_type> frnds_commIDs;
       std::map<graphlab::vertex_id_type, std::vector<graphlab::vertex_id_type> > path_to_DIVas;

       std::vector<single_message> msg2frwd;
       std::vector<std::string> att_group;
       std::vector<size_t> freq;
       

        edge_dir_type gather_edges(icontext_type& context,
                              const vertex_type& vertex) const {
                return graphlab::ALL_EDGES;
        }

        edge_dir_type  scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
                return graphlab::NO_EDGES;
        }

        data_gatherer_directN gather(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
          data_gatherer_directN att_n_v;
          if(context.iteration() ==0){
            if(edge.source().data().join_time <= snapShotID)
               att_n_v.frnds_IDs.push_back(edge.source().id());
          }
          else if(context.iteration() > 0 && context.iteration() < 50 ){
             if(edge.source().data().join_time <= snapShotID){
                att_n_v.frnds_commIDs.insert(std::make_pair(edge.source().id(), edge.source().data().commID));
                att_n_v.path_to_DIVas.insert(std::make_pair(edge.source().id(), edge.source().data().path_to_Diva));
             }
          }
          else if(context.iteration() == 50 ){    
            if(edge.source().data().att_values.size() >0 && edge.source().data().join_time == snapShotID){
               std::multimap<std::string, std::string> att_v;               
               for(std::multimap<std::string, std::string>::iterator iter = edge.source().data().att_values.begin(); 
                   iter!=edge.source().data().att_values.end(); iter++){
                   std::string vk = iter->first;
                   std::string v = iter->second;
                   att_v.insert(std::make_pair(vk, v));
              }
              att_n_v.frnds_profiles.insert(std::make_pair(edge.source().id(),att_v));
           }
          }
          return att_n_v;
        }

        void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
        }


   //====================================================================================

   void init(icontext_type& context,
                const vertex_type& vertex, const data_message &msg) {
        if(context.iteration() >= 1 ){
          for(size_t i=0; i<msg.msgs.size() ; i++){
             single_message reply =  msg.msgs[i];
             //check if I'm central node of message CommID
             if(vertex.data().commID == vertex.id()){
                for(size_t x=0; x< reply.att_group.size(); x++)
                   att_group.push_back(reply.att_group[x]);
                for(size_t x=0; x< reply.freq.size(); x++)
                   freq.push_back(reply.freq[x]);

             }
            //forward message to central node
             else {
               reply.reciever_id = vertex.data().commID; 
               msg2frwd.push_back(reply);           
             }
             }
          }
    
   }


   void apply(icontext_type& context, vertex_type& vertex,const gather_type& total) {
           if(context.iteration() ==0) {
             //Intialize commID and path_to_diva
             vertex.data().commID = vertex.id();
             for(std::vector<graphlab::vertex_id_type>::const_iterator iter = total.frnds_IDs.begin(); iter != total.frnds_IDs.end(); ++iter)
                if(*iter > vertex.data().commID) 
                  vertex.data().commID = *iter;

             vertex.data().path_to_Diva.push_back(vertex.data().commID);
             //if(vertex.data().commID!=vertex.id())
               //vertex.data().path_to_Diva.insert(vertex.data().path_to_Diva.begin(), vertex.id()); 
            }
            else if(context.iteration() >0 && context.iteration()<50){
              graphlab::vertex_id_type connector;

              vertex.data().frnds_commIDs.clear();
              for(std::map<graphlab::vertex_id_type, graphlab::vertex_id_type>::const_iterator iter = total.frnds_commIDs.begin(); iter != total.frnds_commIDs.end(); ++iter)
                if(iter->first != vertex.id())
                   vertex.data().frnds_commIDs.insert(std::make_pair(iter->first, iter->second));

              vertex.data().Divas_paths.clear();
	      for(std::map<graphlab::vertex_id_type, std::vector<graphlab::vertex_id_type> >::const_iterator iter = total.path_to_DIVas.begin(); iter != total.path_to_DIVas.end(); ++iter){
                  vertex.data().Divas_paths.insert(std::make_pair(iter->first, iter->second));
               }
                /*context.cout() << context.iteration() << " \n ";
                for(std::map<graphlab::vertex_id_type, graphlab::vertex_id_type>::iterator iter = vertex.data().frnds_commIDs.begin();
                       iter != vertex.data().frnds_commIDs.end(); ++iter){
                        context.cout() << iter->first << " , " << iter->second<< " \n ";
                }*/

              
               vertex.data().find_DIVas();
                //adjust the path
                if(vertex.data().commID == vertex.id())
                    vertex.data().path_to_Diva.clear();
                else{
                   bool df= false;
                for(std::map<graphlab::vertex_id_type, graphlab::vertex_id_type>::const_iterator iter = total.frnds_commIDs.begin(); iter != total.frnds_commIDs.end(); ++iter)
                if(vertex.data().commID == iter->first){
                   vertex.data().path_to_Diva.clear();
                   vertex.data().path_to_Diva.push_back(iter->first);
                   df = true;
                   break;
                }
                if(!df){
                 //find the shortest path to DIVa
                 int len=-1;
                 std::vector<graphlab::vertex_id_type> path;
                 for(std::map<graphlab::vertex_id_type,std::vector<graphlab::vertex_id_type> >::const_iterator iter = total.path_to_DIVas.begin(); 
                    iter!= total.path_to_DIVas.end(); iter++){
                    if(vertex.data().commID == vertex.data().frnds_commIDs[iter->first]){
                      if(len==-1){
                      connector =  iter->first;
                      path = iter->second;
                      len = path.size();
                      //context.cout() << len <<"\n";
                    }
                    else{
                     size_t l = iter->second.size();
                     if(l<len){
                       connector =  iter->first;
                       path = iter->second;
                       len =l;
                     }
                    }
                 }
                }
                if(count(path, connector) ==0)
                    path.insert(path.begin(), connector);

                vertex.data().path_to_Diva = path;
               }
              }
            }
            else{
             for(std::map< graphlab::vertex_id_type, std::multimap<std::string, std::string> >::const_iterator iter = total.frnds_profiles.begin(); iter!= total.frnds_profiles.end(); ++iter)
               frnds_profiles.insert(std::make_pair(iter->first, iter->second));
               vertex.data().compute_LCAS(frnds_profiles, vertex.data().lCAS, vertex.data().lFreq, vertex.data().lSize); 
               for(std::multimap<graphlab::vertex_id_type, std::vector<std::string> >::iterator it = vertex.data().lCAS.begin(); it != vertex.data().lCAS.end(); it++){
                  single_message msg;
                  msg.sender_id = vertex.id();
                  msg.reciever_id = it->first;
                  msg.path = vertex.data().Divas_paths[it->first];
                  msg.att_group = it->second;
                  msg.freq = vertex.data().lFreq[it->first];
                  msg.size =  vertex.data().lSize[it->first];
                  data_message dmsg;
                  dmsg.msgs.push_back(msg);
                  context.signal_vid(it->first, dmsg);
               }
           }

           context.signal(vertex);
            /*
            //aggregate community results
            if(context.iteration() == 8 ){
               std::vector<size_t> count;
               for(size_t x =0; x< vertex.data().att_group.size(); x++){
                  int index1 = find(vertex.data().att_group_clstr, vertex.data().att_group[x]);
                  if(index1 <0 ){
                     vertex.data().att_group_clstr.push_back(vertex.data().att_group[x]);
                     vertex.data().freq_clstr.push_back(vertex.data().freq[x]);
                     vertex.data().supp_clstr.push_back(vertex.data().CAG_supp[x]);
                     vertex.data().conf_clstr.push_back(vertex.data().CAG_conf[x]);
                  }
                  else{
                     vertex.data().freq_clstr[index1] += vertex.data().freq[x];
                     vertex.data().supp_clstr[index1] += vertex.data().CAG_supp[x];
                     vertex.data().conf_clstr[index1] += vertex.data().CAG_conf[x];
                  }
               }
            }
            else if(context.iteration() > 0 ){
              //check for messages to forward
             if(msg2frwd.size()>0)
               for(size_t i=0; i<msg2frwd.size(); i++){
                  data_message msg;
                  msg.msgs.push_back(msg2frwd[i]);
                  context.signal_vid(vertex.data().commID, msg);
               }
                 //aggregate with the new recieved msgs
                for(size_t x =0; x< att_group.size(); x++){
                      vertex.data().att_group.push_back(att_group[x]);
                      vertex.data().freq.push_back(freq[x]);
                      vertex.data().CAG_supp.push_back(CAG_supp[x]);
                      vertex.data().CAG_conf.push_back(CAG_conf[x]);
                   }

               //signal myself to keep alive to write done my results @ the end
               context.signal_vid(vertex.id());
           }
           //intialization
           else {
           //first loop over gathered data and add them to vertex att_values
            //context.cout() << "in apply \n";
            //claculating similarity
            std::vector<std::string> all_values;
            std::vector< std::multimap<std::string, std::string> > frnds_profiles;
            std::vector<std::string> att_names;

           att_names.push_back("school");
           att_names.push_back("major");
           att_names.push_back("employer");
           att_names.push_back("places_lived");

            //
            for(std::vector< std::multimap<std::string, std::string> >::const_iterator iter = total.frnds_profiles.begin(); iter!= total.frnds_profiles.end(); ++iter)
               frnds_profiles.push_back(*iter);

           //if( vertex.data().att_values.size()>0)
           {
           //calculting the frequencies using node's att_values
           for(size_t i=0;i<att_names.size()-1; i++){
               std::string a_n1 = att_names[i];
               std::vector<std::string> triple;
               for(size_t j=i+1; j<att_names.size(); j++){
                  std::string a_n2 = att_names[j];   
                  if(a_n1 != a_n2){
                   for(size_t c=0; c<frnds_profiles.size(); c++){
                   std::multimap<std::string, std::string> frd_p = frnds_profiles[c];
                  //context.cout() <<"an1: "<<a_n1 << " with an2: " << a_n2 << "\n";
                  std::pair<std::multimap<std::string, std::string>::iterator,std::multimap<std::string, std::string>::iterator > 
                       iter_n1 = frd_p.equal_range(a_n1);

                  std::pair<std::multimap<std::string, std::string>::iterator,std::multimap<std::string, std::string>::iterator >
                       iter_n2 = frd_p.equal_range(a_n2);

                  std::vector<std::string> val_n1, val_n2, f_val_n1, f_val_n2;
                  for (std::multimap<std::string, std::string>::iterator it=iter_n1.first; it!=iter_n1.second; ++it){
                      val_n1.push_back(it->second);
                      //context.cout() <<"values for an1: " << it->second << "\n";
                  }
                  for (std::multimap<std::string, std::string>::iterator it=iter_n2.first; it!=iter_n2.second; ++it){
                      val_n2.push_back(it->second);
                      //context.cout() <<"values for an2: " << it->second << "\n";
                  }
                  std::string value1, value2;
                  size_t len1 = val_n1.size();
                  size_t len2 = val_n2.size();
                  std::vector<double> val_supp;
                  for(size_t x = 0; x<len1 ;x++){
                     size_t count =1;
                     for(size_t y = x+1; y<len1 ;y++){
                        if(val_n1[x] == val_n1[y])
                          count ++;
                     }
                     double supp = double(count)/double(len1+len2);
                     //supp /=100;
                     val_supp.push_back(supp);
                  }
                  for(size_t x = 0; x<len1 ;x++){
                       value1 = val_n1[x];
                       size_t freq_n=1;
                       //double sup_unit=double(1)/double(len1+len2);
                       //sup_unit /=100;
                       for(size_t y = 0; y<len2 ;y++){
                          value2 = val_n2[y];
                           std::string key1 = a_n1 +  "&" + a_n2;

               for(size_t rest_of_frnds = 0;rest_of_frnds<frnds_profiles.size(); rest_of_frnds ++){
                if(rest_of_frnds != c){
                  std::pair<std::multimap<std::string, std::string>::iterator,std::multimap<std::string, std::string>::iterator >
                       iter_n1 = frnds_profiles[rest_of_frnds].equal_range(a_n1);
                  std::pair<std::multimap<std::string, std::string>::iterator,std::multimap<std::string, std::string>::iterator >
                       iter_n2 = frnds_profiles[rest_of_frnds].equal_range(a_n2);
                  std::vector<std::string> f_val_n1, f_val_n2;
                  for (std::multimap<std::string, std::string>::iterator it=iter_n1.first; it!=iter_n1.second; ++it){
                      f_val_n1.push_back(it->second);
                  }
                  for (std::multimap<std::string, std::string>::iterator it=iter_n2.first; it!=iter_n2.second; ++it){
                      f_val_n2.push_back(it->second);
                  }
                
                  if(find(f_val_n1, value1) >0 && find(f_val_n2, value2) >0){
                      for(size_t inf1=0; inf1< f_val_n1.size(); inf1++)
                               if(value1 == f_val_n1[inf1])
                                   freq_n++;
                           for(size_t inf1=0; inf1< f_val_n2.size(); inf1++)
                               if(value2 == f_val_n2[inf1])
                                   freq_n++;

                  }   
                }
               } 
                            

                           int index1 = find(vertex.data().att_group, key1);

                           if(index1 <0 ){
                             vertex.data().att_group.push_back(key1);
                             vertex.data().freq.push_back(freq_n);
                             vertex.data().CAG_supp.push_back(len1+len2);
                             vertex.data().CAG_conf.push_back(val_supp[x]);
                           }
                           else{
                                  vertex.data().freq[index1] ++;
                                  //vertex.data().CAG_supp[index1] +=sup_unit;
                                  //vertex.data().CAG_conf[index1] = vertex.data().CAG_supp[index1]/val_supp[x];
                         } 
                       }
                    }  
                  }
                }
               }
              }
              single_message reply;
              reply.sender_id = vertex.id();
              reply.reciever_id =  vertex.data().commID;
              reply.att_group = vertex.data().att_group;
              reply.freq = vertex.data().freq;
              reply.supp =  vertex.data().CAG_supp;
              reply.conf = vertex.data().CAG_conf;
              reply.size.push_back( vertex.data().friends.size());
              data_message msg;
              msg.msgs.push_back(reply);
              context.signal_vid(vertex.data().commID, msg);
           }
        }*/
 }
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//############################ line parser to load data from file ##################################################//

bool line_parser(graph_type& graph,
                   const std::string& filename,
                   const std::string& textline)
 {
    if(textline.empty()) return true;

    std::stringstream strm(textline);

    graphlab::vertex_id_type vid;
    std::vector<graphlab::vertex_id_type> frnds, c_frnds;
    std::string attribute, token, token_sub, token_list;
    bool  hasFrnds = false, hasAttrib=false;
    size_t join_time;
    std::vector<std::string> att_names;
    std::multimap<std::string, std::string>  att_values,c_att_values;

    att_names.push_back("school");
    att_names.push_back("major");
    att_names.push_back("employer");
    att_names.push_back("places_lived");

/////////////////////////////////////////////////Load data    
    strm >> vid;
    if(vid!=0){
      graph.add_vertex(vid, vdata(join_time, att_values));
    strm >> attribute;
    std::istringstream ss8(attribute); 
    while(std::getline(ss8, token, '#')){
       size_t fID = atoi(token.c_str());
       if(fID !=0){
       graph.add_edge(fID, vid);
       graph.add_edge(vid, fID);}
    }
   }

   /*
    strm >> vid;
    strm >> join_time;
    
    strm >> attribute;
    if(attribute != "null"){
       hasFrnds = true;
       std::istringstream ss(attribute);
       while(std::getline(ss, token, '\\')) {
          size_t fID = atoi(token.c_str());
          if(fID != vid){
             frnds.push_back(fID);
          }
       }
    }
   
   strm >> attribute;
  if(attribute != "null"){ 
  hasAttrib = true;
  std::istringstream ss8(attribute);

   while(std::getline(ss8, token, '@')){
      int index =0;
      std::string type, values;
      std::istringstream ss10(token);
      while(std::getline(ss10, token_sub, '#')){
         if(index ==0){
            values = token_sub;
            index ++;
         }
         else{
             type = token_sub;
             std::istringstream ss11(values);
             while(std::getline(ss11, token_list, '\\')){   
                 if(token_list!= "")
	            att_values.insert(std::make_pair(type,token_list));  
             }
         }
      }
   }
   }

///////////////////////////////////////////Add node to the graph
  if(hasAttrib || hasFrnds) {  
   graph.add_vertex(vid, vdata(join_time, att_values));
    //if(isCent)
      for(int  i=0; i< frnds.size(); i++)
         graph.add_edge(frnds[i], vid);
   }
   */

   return true;
 }


//#################################  Graph Writer ######################################//
//#################################  Graph Writer ######################################//
struct graph_writer {
  std::string save_vertex(graph_type::vertex_type v) {


        std::stringstream strm;
        std::vector<std::string> att_vec;
        std::vector<size_t> att_freq;
        strm << "============================================== \n";
        strm <<  v.id() << " , " << v.data().commID << " \n ";
        for(std::map<graphlab::vertex_id_type, size_t>::iterator iter = v.data().clusters.begin();
                       iter != v.data().clusters.end(); ++iter){
                        strm << iter->first << " , " << iter->second<< " \n ";
                }        

        strm << " \n "; 
        for(size_t i=0; i<v.data().path_to_Diva.size(); i++)
            strm <<v.data().path_to_Diva[i] <<" , " ; 
        strm << "\n";

        for(std::map<graphlab::vertex_id_type, std::vector<graphlab::vertex_id_type> >::iterator iter = v.data().Divas_paths.begin();
               iter != v.data().Divas_paths.end(); iter++){
                  strm << iter->first << "( ";
                  
                  for(size_t i=0; i< iter->second.size(); i++)
                      strm << iter->second[i] << " , ";
                 strm << ")\n";
            }

        /*if(v.data().att_group_clstr.size()>0){    
          strm << "============================================== \n";
          strm <<  v.id() << " , " << v.data().friends.size() << "\n";
          
          for(size_t i=0; i< v.data().att_group_clstr.size(); i++){
             double supp = double(v.data().freq_clstr[i])/double(v.data().supp_clstr[i]);
             double conf = supp/double(v.data().conf_clstr[i]);
             strm << v.data().att_group_clstr[i] << " , " << v.data().freq_clstr[i] << " , " << supp << " , " << conf << "\n";
          }
       }*/
       //for(size_t i=0; i< v.data().att_group.size(); i++)
            // strm << v.data().att_group[i] << " , " << v.data().freq[i] << "\n";

    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
};

//######################################## main function #########################################//
//######################################## main function #########################################//
  int main(int argc, char** argv) {

        graphlab::command_line_options clopts("Dosn: fill neighbors map");
        std::string file_name;
        clopts.attach_option("graph", file_name, "Graph input. reads all graphs matching prefix*");

        if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
        if (file_name == "") {
                std::cout << "--graph is not optional\n";
                clopts.print_description();
                return EXIT_FAILURE;
        }

        // Initialize control plane using mpi
        graphlab::mpi_tools::init(argc, argv);
        graphlab::distributed_control dc;
        graphlab::launch_metric_server();

        // load graph
        graph_type graph(dc, clopts);      //graph_type graph(dc);
         
        graph.load(file_name, line_parser); //load_from_hdfs("hdfs://graphLabLarge1:154310/data/graph.txt", line_parser);
        //graph.load_format(file_name, "tsv");
        graph.finalize();
        dc.cout() << "Number of vertices: " << graph.num_vertices() << std::endl
              << "Number of edges:    " << graph.num_edges() << std::endl;

   
   size_t MAX_ITERATIONS=3;
   clopts.get_engine_args().set_option("max_iterations", MAX_ITERATIONS);
   graphlab::omni_engine<DIVa_NodeCentric> engine2(dc, graph, "sync", clopts);
   engine2.signal_all();
   engine2.start();
   std::cout << "Runtime: " << engine2.elapsed_seconds();
   graph.save(file_name, graph_writer(),
             false,    // do not gzip
             true,     // save vertices
             false);   // do not save edges
  
   graphlab::mpi_tools::finalize();



   return EXIT_SUCCESS;
  }



