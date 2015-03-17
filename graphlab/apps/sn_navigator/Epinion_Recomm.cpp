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
//////////////////////////////////////////////////////////
//gossip message exchanged between social graph nodes
//////////////////////////////////////////////////////////

struct single_message {
    graphlab::vertex_id_type sender_id;
    graphlab::vertex_id_type reciever_id;  
    std::multimap<graphlab::vertex_id_type, double>  att_values;
    std::vector<std::string> att_group;
    std::vector<size_t> freq;

    single_message(): att_values() { }

   void save(graphlab::oarchive& oarc) const {
                oarc <<sender_id<<reciever_id;
                size_t num = att_values.size();
                oarc << num;
                for(std::multimap<graphlab::vertex_id_type, double>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                }
               num = att_group.size();
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


   }
   void load(graphlab::iarchive& iarc) {
          iarc>>sender_id>>reciever_id;
          size_t num,f;
          graphlab::vertex_id_type vk;
          std::string s;
          double v;
          iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
               }
         iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> s;
                  att_group.push_back(s);
               }
        iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> f;
                  freq.push_back(f);
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
        bool is_centroid;        
         graphlab::vertex_id_type vid;
        std::vector<graphlab::vertex_id_type> friends;
        std::vector<std::string> att_names;
        std::vector<double> products;
        std::multimap<graphlab::vertex_id_type, double>  att_values;
        std::multimap<std::string, size_t > att_v_n_freq;
        std::vector<std::string> att_group;
        std::vector<size_t> freq;
        std::vector<std::string> att_group_clstr;
        std::vector<size_t> freq_clstr;

     vdata() :
                friends(), att_names(), att_values(), att_v_n_freq(), att_group(),freq(), att_group_clstr(),freq_clstr(), is_centroid(false) {}

      explicit vdata(graphlab::vertex_id_type id ,std::vector<graphlab::vertex_id_type> f,std::vector< std::string> a_n, std::multimap<graphlab::vertex_id_type, double> a_v,  bool isc, std::vector< double> p) :
               vid (id), friends(f),att_names(a_n), att_values(a_v), att_v_n_freq(), att_group(),freq(), att_group_clstr(),freq_clstr(), is_centroid(isc), products(p) {}

        void save(graphlab::oarchive& oarc) const {
               
                //store list of friends
                oarc << vid;
                oarc << is_centroid;

                size_t num = friends.size();
                oarc << num;
                for(size_t i = 0;i < num; ++i){
                        oarc << friends[i];
                }

                //CbIV
                num = att_names.size();
                oarc << num;
                for(std::vector<std::string>::const_iterator iter = att_names.begin();
                       iter != att_names.end(); ++iter){
                               oarc << *iter;
                }
                 num = products.size();
                oarc << num;
                for(size_t in=0; in<num; in ++){
                     //if(is_centroid && in <5)
                     //std::cout << "out vid:" << vid <<  ",  " <<products[in] << "\n";
                     oarc  << products[in];
                }
                /*for(std::vector<std::string>::const_iterator iter = products.begin();
                       iter != products.end(); ++iter){
                               oarc << *iter;
                }*/

                num = att_values.size();
                oarc << num;
                for(std::multimap<graphlab::vertex_id_type, double>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                       //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
                }
                num = att_group.size();
                oarc << num;
                for(std::vector<std::string>::const_iterator iter = att_group.begin();
                       iter != att_group.end(); ++iter){
                               oarc << *iter;
                }
                num = freq.size();
                oarc << num;
                for(std::vector<std::size_t>::const_iterator iter = freq.begin();
                       iter != freq.end(); ++iter){
                               oarc << *iter;
                }
                num = att_group_clstr.size();
                oarc << num;
                for(std::vector<std::string>::const_iterator iter = att_group_clstr.begin();
                       iter != att_group_clstr.end(); ++iter){
                               oarc << *iter;
                }
                num = freq_clstr.size();
                oarc << num;
                for(std::vector<std::size_t>::const_iterator iter = freq_clstr.begin();
                       iter != freq_clstr.end(); ++iter){
                               oarc << *iter;
                }
                num = att_v_n_freq.size();
                oarc << num;
                for(std::multimap<std::string, size_t>::const_iterator iter = att_v_n_freq.begin();
                       iter != att_v_n_freq.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                               /*std::vector<std::string> vlist =  iter->second;
                               for(size_t vi =0; vi<vlist.size(); vi++)
                                   oarc << vlist[vi];*/
                }
                



        }
     void load(graphlab::iarchive& iarc) {
               
               std::string v,vk;
               graphlab::vertex_id_type gl_id;
               size_t num = 0;
               double val;
               iarc >> vid;
               iarc >> is_centroid;

               iarc >> num;
                for(size_t i = 0; i < num; i++) {
                  graphlab::vertex_id_type f_id = 0;
                  iarc >> f_id;
                  friends.push_back(f_id);
                }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk ;
                  att_names.push_back(vk);
               }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> val ;
                  products.push_back(val);
                  //if(is_centroid && i <5)
                     //std::cout << "in: vid: " << vid << " , " <<products[i] << "\n";

               }


               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> gl_id >> val;
                  att_values.insert(std::make_pair(gl_id,val));
                  //if(vk != "interest")
                  //std::cout <<"key: "<<vk<<", val: "<<v<< "\n";
               }

               iarc >> num;
               size_t f=0;
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
                  iarc >> vk;
                  att_group_clstr.push_back(vk);
               }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> f;
                  freq_clstr.push_back(f);
               }
               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk ;
                  iarc >> f ;
                  //Att1value
                  /*std::vector<std::string> vs;
                  iarc >> v;
                  vs.push_back(v);
                  //Att2value
                  iarc >> v;
                  vs.push_back(v);
                  //Att3value
                  iarc >> v;
                  vs.push_back(v);*/

                  att_v_n_freq.insert(std::make_pair(vk,f));
               }
        }

 };


///////////////////////////////global functions/////////////////////////////////

int find(std::vector<double> vec, double key){
    int index =-1;
    
    for(size_t i=0; i<vec.size(); i++){
       if(vec[i] == key){
          index =i; 
          break;
       }
    }
    return index;
}

int findS(std::vector<std::string> vec, std::string key){
    int index =-1;
    size_t len=0;
    double key_val1, key_val2, list_val1, list_val2;
    std::string token;
   
    for(size_t i=0; i<vec.size(); i++){
         if(vec[i] == key){
          index =i;
          break;
       }
    }
    return index;
}


std::vector<graphlab::vertex_id_type> getRandomCent( graphlab::vertex_id_type vid ){
   std::vector<graphlab::vertex_id_type> val;

   std::multimap<graphlab::vertex_id_type, graphlab::vertex_id_type> centroids;
   //164//
   centroids.insert(std::make_pair(164, 16));
   centroids.insert(std::make_pair(164,124 ));
   //294//
   centroids.insert(std::make_pair(294,0 ));
   centroids.insert(std::make_pair(294, 16));
   //393//
   centroids.insert(std::make_pair(393, 0));
   centroids.insert(std::make_pair(393, 16));
   //0//
   centroids.insert(std::make_pair(0, 294));
   centroids.insert(std::make_pair(0, 393));
   //16//
   centroids.insert(std::make_pair(16, 124));
   centroids.insert(std::make_pair(16, 164));
   //215//
   centroids.insert(std::make_pair(215, 654));
   centroids.insert(std::make_pair(215, 0));
   //258//
   centroids.insert(std::make_pair(258, 124));
   centroids.insert(std::make_pair(258, 164));
   //124//
   centroids.insert(std::make_pair(124, 16));
   centroids.insert(std::make_pair(124, 164));
   //654//
   centroids.insert(std::make_pair(654, 0));
   centroids.insert(std::make_pair(654, 215));

    std::pair<std::multimap<graphlab::vertex_id_type, graphlab::vertex_id_type>::iterator,std::multimap<graphlab::vertex_id_type, graphlab::vertex_id_type>::iterator >
                       iter = centroids.equal_range(vid);
                  
                  for (std::multimap<graphlab::vertex_id_type, graphlab::vertex_id_type>::iterator it=iter.first; it!=iter.second; ++it){
                      val.push_back(it->second);
                      std::cout<< it->second <<"\n";
                  }

   return val;
   
}

size_t getIntraConnCount ( graphlab::vertex_id_type vid ){
    size_t links = 0;

   std::map< graphlab::vertex_id_type,size_t> centroids;
   std::map< graphlab::vertex_id_type,size_t>::iterator it;
   centroids[164] = 4;
   centroids[294] = 2;
   centroids[393] = 2;
   centroids[0] = 5;
   centroids[16] = 3;
   centroids[215] = 7;
   centroids[258] = 3;
   centroids[124] = 3;
   centroids[654] = 1;
   
   it=centroids.find(vid);
   links  = it->second; 
     

    return links;
}

bool isCentroid(  graphlab::vertex_id_type vid ){
   bool chk = false;
   std::map< graphlab::vertex_id_type,  graphlab::vertex_id_type> centroids;
   std::map< graphlab::vertex_id_type,  graphlab::vertex_id_type>::iterator it;
   //0, 16, 124, 164, 215, 258, 294, 393, and 654
   centroids[8319] = 8319;
   centroids[13677] = 13677;
   centroids[5399] = 5399;
   centroids[18220] = 18220;
   centroids[20240] = 20240;
   centroids[9831] = 9831;
   centroids[6876] = 6876;
   centroids[2760] = 2760;
   centroids[5369] = 5369;

   it=centroids.find(vid);
   if(it != centroids.end())
     chk = true;

   return chk;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//distributed graph definition
typedef graphlab::distributed_graph<vdata, edata> graph_type;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct empty_data_gatherer {
        empty_data_gatherer (){}
        void save(graphlab::oarchive& oarc) const {

        }
        void load(graphlab::iarchive& iarc) {

        }
       empty_data_gatherer& operator+=(const empty_data_gatherer& other) {
              return *this;
         }
};

class collect_attributes:public graphlab::ivertex_program<graph_type,empty_data_gatherer,data_message>, public graphlab::IS_POD_TYPE {

public:
 
       std::multimap<graphlab::vertex_id_type, double> att_values;
       std::vector<single_message> msg2frwd;
       std::vector<std::string> att_group_clstr;
       std::vector<size_t> freq_clstr;

        edge_dir_type gather_edges(icontext_type& context,
                              const vertex_type& vertex) const {
                return graphlab::ALL_EDGES;
        }

        edge_dir_type  scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
                return graphlab::NO_EDGES;
        }

        empty_data_gatherer gather(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
               empty_data_gatherer edg;
               return edg;
        }

        void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
        }


   //====================================================================================

   void init(icontext_type& context,
                const vertex_type& vertex, const data_message &msg) {

       if(context.iteration() == 1 ){
          for(size_t i=0; i<msg.msgs.size() ; i++){
             if(vertex.id() == msg.msgs[i].reciever_id){
                single_message reply;
                reply.sender_id = vertex.id();
                reply.reciever_id = msg.msgs[i].sender_id;
                 for(std::vector< double>::const_iterator iter = vertex.data().products.begin();
                   iter!=vertex.data().products.end(); iter++){
                   double v = *iter;
                   reply.att_values.insert(std::make_pair(vertex.id(), v));
                  }
                msg2frwd.push_back(reply);
             }
          }        
          
       }
       else if(context.iteration() == 2 ){
          for(size_t i=0; i<msg.msgs.size() ; i++){
             if(vertex.id() == msg.msgs[i].reciever_id){
                 single_message reply =  msg.msgs[i];
                 for(std::multimap<graphlab::vertex_id_type, double>::iterator iter = reply.att_values.begin();
                   iter!=reply.att_values.end(); iter++){
                   graphlab::vertex_id_type vk = iter->first;
                   double v = iter->second;
                   att_values.insert(std::make_pair(vk, v));
              }

             }
          }
       }
       else if(context.iteration() == 3 ){
            msg2frwd.clear();
            for(size_t i=0; i<msg.msgs.size() ; i++){
             if(vertex.id() == msg.msgs[i].reciever_id){
                single_message reply;
                reply.sender_id = vertex.id();
                reply.reciever_id = msg.msgs[i].sender_id;
                reply.att_group = vertex.data().att_group;
                reply.freq = vertex.data().freq;
                msg2frwd.push_back(reply);
             }
          }

       }
     else if(context.iteration() == 4 ){
         for(size_t i=0; i<msg.msgs.size() ; i++){
             if(vertex.id() == msg.msgs[i].reciever_id){
                 single_message reply =  msg.msgs[i];
                 for(std::vector<std::string>::iterator iter = reply.att_group.begin();
                   iter!=reply.att_group.end(); iter++){
                   att_group_clstr.push_back(*iter);
                   //std::cout <<"CAG: "<<*iter<< "\n";
                 }
                for(std::vector<size_t>::iterator iter = reply.freq.begin();
                   iter!=reply.freq.end(); iter++){
                   freq_clstr.push_back(*iter);
                   //std::cout <<"freq: "<<*iter<< "\n";
                 }
             }
          }

     }
    
   }


   void apply(icontext_type& context, vertex_type& vertex,const gather_type& total) {
         
        //context.cout() << "apply it: " << context.iteration() << "\n" ;
        if(context.iteration() == 0){
           if(isCentroid( vertex.id())){
              for(std::vector<graphlab::vertex_id_type>::iterator iter = vertex.data().friends.begin(); iter!= vertex.data().friends.end(); iter++){
                  single_message req;
                  req.sender_id = vertex.id();
                  req.reciever_id = *iter;
                  data_message msg;
                  msg.msgs.push_back(req);
                  context.signal_vid(*iter, msg);
              }
           }
           context.signal_vid(vertex.id()); 
        }
        else if(context.iteration() == 1){
           for(size_t i=0; i<msg2frwd.size(); i++){
                  data_message msg;
                  msg.msgs.push_back(msg2frwd[i]);
                  context.signal_vid(msg2frwd[i].reciever_id, msg);
           }
        }

        else if(context.iteration() == 2 && isCentroid( vertex.id()) ){
           //calculting the frequencies using node's att_values
           std::vector<double> val_n1, val_n2;
           std::stringstream convert, c1;
           //get node's products first
           //context.cout() << "node id " << vertex.id() << " ,prosuct  set size" <<  vertex.data().products.size() << "\n";
           for(size_t in=0; in<  vertex.data().products.size(); in++) {
             if(find (val_n1, vertex.data().products[in] )==-1) {
                 val_n1.push_back( vertex.data().products[in]);
                 //context.cout() << "n" << vertex.id() << " , e" << in <<" : " <<  vertex.data().products[in] << "\n";
              }

            }

           //insert direc friends products
           for(std::multimap<graphlab::vertex_id_type, double>::iterator iter = att_values.begin(); iter!= att_values.end(); ++iter){
               vertex.data().att_values.insert(std::make_pair(iter->first, iter->second));
               //if(vertex.id() == 16)
               //context.cout() << "adding " << iter->first <<" , " <<  iter->second << "\n";
               //context.cout() << "n" << vertex.id() << " , e" << in <<" : " <<   iter->second << " , f" << iter->first << "\n";
           }
           context.cout() << "\n=======================================================\n ";
           
           for(size_t i=0; i<vertex.data().friends.size(); i++){
               graphlab::vertex_id_type id = vertex.data().friends[i];
               //if(vertex.id() == 16)
                   //context.cout() << "checking with friend " <<  id << "\n";
               std::pair<std::multimap<graphlab::vertex_id_type, double>::iterator,std::multimap<graphlab::vertex_id_type, double>::iterator >
                   iter_f_prod = vertex.data().att_values.equal_range(id);
               for (std::multimap<graphlab::vertex_id_type, double>::iterator it=iter_f_prod.first; it!=iter_f_prod.second; ++it) {
                   if(find(val_n2, it->second) ==-1){
                      val_n2.push_back(it->second);
                      //if(vertex.id() == 16)
                         //context.cout() << " e:"  <<   it->second << " , f" << it->first << "\n";
                   }
               }
               //context.cout() << "vid: " << vertex.id() << "1st set size" <<  val_n1.size() << ",2nd set size" <<  val_n2.size() << "\n";
               double value1, value2;
               size_t len = val_n1.size();
                               
                for(size_t p1=0; p1<len; p1++){
                    for(size_t p2=p1+1; p2<len; p2++){
                       value1 = val_n1[p1];
                       value2 = val_n1[p2];
                       //context.cout() << "comparing " << value1 <<" , " <<  value2 << "\n";
                       if(value1 != value2)
                       if((find(val_n2, value1) >=0) && (find(val_n2, value2) >=0)){
                              //context.cout() << "adding " << value1 <<" , " <<  value2 << "\n";
                              std::stringstream c1,c2;
                              std::string key;
                              c1 << value1;
                              c2 << value2;
                              key = "c" + c1.str() + "_#_c" + c2.str();
                              int index = findS(vertex.data().att_group, key);
                              if(index ==-1){
                                  vertex.data().att_group.push_back(key);
                                  vertex.data().freq.push_back(1);
                               }
                               else {
                                   //context.cout() <<"vid : " << vertex.id() << " , index of key " << index << " , key: "<< key<< "\n";
                                   vertex.data().freq[index] += 1;
                               }
                                  
                       }                 
                    }
                }
            } 
                
           context.cout() << vertex.id() << " , " << vertex.data().att_group.size() << "\n";
            //context.cout() << "************************ node finished " <<vertex.id()  << "\n";
            /*for(size_t index=0; index < vertex.data().freq.size(); index ++){
                vertex.data().freq[index] = (double(vertex.data().freq[index]) / double(vertex.data().friends.size()))*100;
                 vertex.data().freq[index] *= getIntraConnCount(vertex.id());
            }*/
           //pick 3 random nodes and forward discovered CAG to them.
          std::vector<graphlab::vertex_id_type> rndcentids = getRandomCent(vertex.id());
          for(size_t send=0; send<rndcentids.size() ; send ++) {
             graphlab::vertex_id_type rndCent = rndcentids[send];
             single_message req;
             req.sender_id = vertex.id();
             req.reciever_id = rndCent;
             data_message msg;
             msg.msgs.push_back(req);
             context.signal_vid(rndCent, msg);
          }
          context.signal_vid(vertex.id());
         
        }
       else if(context.iteration() == 3 ){
            for(size_t i=0; i<msg2frwd.size(); i++){
                  data_message msg;
                  msg.msgs.push_back(msg2frwd[i]);
                  context.signal_vid(msg2frwd[i].reciever_id, msg);
           }
           context.signal_vid(vertex.id());
       }
     else if(context.iteration() == 4 ){
         //sum aggregated CAGs
         size_t length = att_group_clstr.size();
         if(freq_clstr.size() < length)
            length = freq_clstr.size();
         for(size_t index =0; index < length; index++){
            std::string key = att_group_clstr[index];
            int keyindex = findS(vertex.data().att_group_clstr, key);
            if(keyindex==-1){
               vertex.data().att_group_clstr.push_back(key);
               vertex.data().freq_clstr.push_back(freq_clstr[index]);
            }
            else
              vertex.data().freq_clstr[keyindex] += freq_clstr[index];
            
         }
         context.cout() << "************************ node finished " <<vertex.id()  << "\n";
     }
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
    std::vector<graphlab::vertex_id_type> frnds;
    std::string attribute, token, token_sub, token_list;
    bool isCent = false;
    std::vector<std::string> att_names;
    std::vector<double> num_values;
    std::multimap<graphlab::vertex_id_type, double>  att_values;

    //att_names.push_back("vid");
    //att_names.push_back("category");
    //att_names.push_back("rating");

/////////////////////////////////////////////////Load data    
    strm >> vid;
    strm >> attribute;
   
    //graph.add_vertex(vid);
    isCent = isCentroid (vid);
    if(isCent) {
       std::istringstream ss(attribute);
       while(std::getline(ss, token, '#')) {
          size_t fID = atoi(token.c_str());
          if(fID != vid){
             //std::cout << "node: " << vid  << " , friend: " << fID<< "\n";
             frnds.push_back(fID);
             //graph.add_edge(fID, vid);
          }
       }
    }
  //products
   strm >> attribute;
   //categories
   strm >> attribute;
   std::istringstream ss1(attribute);
   
   while(std::getline(ss1, token, '#')){
      double num_val = atof(token.c_str());
      if(find (num_values, num_val) ==-1){
        
       num_values.push_back(num_val);
       //att_values.insert(std::make_pair(vid,num_val));
     }
   }
 

///////////////////////////////////////////Add node to the graph
    
    graph.add_vertex(vid, vdata(vid,frnds,att_names,att_values,isCent, num_values));
    if(isCent){
       //std::cout << vid << ", p count: " << att_values.size() << ", & values size," << values.size()  <<  "\n"; 
      for(int  i=0; i< frnds.size(); i++)
         graph.add_edge(frnds[i], vid);
    } 
   

   return true;
 }


//#################################  Graph Writer ######################################//
//#################################  Graph Writer ######################################//
struct graph_writer {
  std::string save_vertex(graph_type::vertex_type v) {


        std::stringstream strm;
        
        if(isCentroid(v.id())){
         strm << "============================================== \n";
         strm <<  v.id() << "\n";
         //for(size_t i=0; i< v.data().products.size(); i++)
            //strm<<"e" << i <<" : " <<  v.data().products[i] << "\n";
         /*size_t i=0;
         for(std::multimap<graphlab::vertex_id_type, double>::iterator iter = v.data().att_values.begin(); iter!= v.data().att_values.end() && v.id()==16; ++iter){
                strm<<"e" << i <<" : " <<  iter->second << ", f" << iter->first << "\n";
                i++;
         }*/
         size_t len = v.data().att_group.size();
         if(v.data().freq.size() <len)
             len =v.data().freq.size();
          strm << len << "\n";
         for(size_t i=0; i< len; i++)
             //strm <<  v.data().att_group_clstr[i] <<"," <<  v.data().freq_clstr[i] << "\n";
             strm <<  v.data().att_group[i] <<"," <<  v.data().freq[i] << "\n";
        }
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
   graphlab::omni_engine<collect_attributes> engine2(dc, graph, "sync", clopts);
   engine2.signal_all();
   engine2.start();
   graph.save("EpRec_Data.txt", graph_writer(),
             false,    // do not gzip
             true,     // save vertices
           false);   // do not save edges

   graphlab::mpi_tools::finalize();



   return EXIT_SUCCESS;
}
