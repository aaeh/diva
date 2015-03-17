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

size_t VIEW_SIZE = 50;
size_t PAYLOAD_SIZE = 10;


//////////////////////////////////////////////////////////
//gossip message exchanged between social graph nodes
//////////////////////////////////////////////////////////

struct single_message {
    graphlab::vertex_id_type sender_id;
    graphlab::vertex_id_type reciever_id;  
    std::multimap<std::string, std::string>  att_values;
    std::vector<std::string> att_group;
    std::vector<size_t> freq;
    std::vector<size_t> size;
    std::vector<double> supp;
    std::vector<double> conf;

    single_message(): att_values(),att_group(),freq(),size(), supp(), conf() { }

   void save(graphlab::oarchive& oarc) const {
                oarc <<sender_id<<reciever_id;
                size_t num = att_values.size();
                oarc << num;
                for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin();
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
                num = size.size();
                oarc << num;
                for(std::vector<size_t>::const_iterator iter = size.begin();
                       iter != size.end(); ++iter){
                               oarc << *iter;
                }

                num = supp.size();
                oarc << num;
                for(std::vector<double>::const_iterator iter = supp.begin();
                       iter != supp.end(); ++iter){
                               oarc << *iter;
                }
  
                num = conf.size();
                oarc << num;
                for(std::vector<double>::const_iterator iter = conf.begin();
                       iter != conf.end(); ++iter){
                               oarc << *iter;
                }



   }
   void load(graphlab::iarchive& iarc) {
          iarc>>sender_id>>reciever_id;
          size_t num,f;
          std::string vk,v;
          double val;
          iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
               }

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
                  supp.push_back(f);
               }
         iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> val;
                  conf.push_back(f);
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
                
        size_t commID;
        std::vector<graphlab::vertex_id_type> friends;
        std::vector<std::string> att_names;
        std::multimap<std::string, std::string>  att_values;
        std::vector<std::string> att_group;
        std::vector<size_t> freq;
        std::vector<double> CAG_supp;
        std::vector<double> CAG_conf;
        std::vector<std::string> att_group_clstr;
        std::vector<size_t> freq_clstr;
        std::vector<size_t> size_clstr;
        std::vector<double> supp_clstr;
        std::vector<double> conf_clstr;
        

     vdata() : 
                friends(), att_names(), att_values(), att_group(),freq(), att_group_clstr(),freq_clstr(), CAG_supp(),CAG_conf(), commID(0),size_clstr() {}

      explicit vdata(std::vector<graphlab::vertex_id_type> f,std::vector< std::string> a_n, std::multimap<std::string, std::string> a_v,  size_t cid) :
          friends(f),att_names(a_n), att_values(a_v), att_group(),freq(), att_group_clstr(),freq_clstr() ,CAG_supp(), CAG_conf(), commID(cid), size_clstr() {}

        void save(graphlab::oarchive& oarc) const {
              
                //store list of friends
                
                oarc << commID;
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
                
                num = att_values.size();
                oarc << num;
                for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin();
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
                num = CAG_supp.size();
                oarc << num;
                for(std::vector<double>::const_iterator iter = CAG_supp.begin(); iter != CAG_supp.end(); ++iter){
                               oarc << *iter;
                }
                num = CAG_conf.size();
                oarc << num;
                for(std::vector<double>::const_iterator iter = CAG_conf.begin(); iter != CAG_conf.end(); ++iter){
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
                for(std::vector<size_t>::const_iterator iter = freq_clstr.begin();
                       iter != freq_clstr.end(); ++iter){
                               oarc << *iter;
                }
                num = size_clstr.size();
                oarc << num;
                for(std::vector<std::size_t>::const_iterator iter = size_clstr.begin();
                       iter != size_clstr.end(); ++iter){
                               oarc << *iter;
                }
                num = supp_clstr.size();
                oarc << num;
                for(std::vector<double>::const_iterator iter = supp_clstr.begin();
                       iter != supp_clstr.end(); ++iter){
                               oarc << *iter;
                }
                num = conf_clstr.size();
                oarc << num;
                for(std::vector<double>::const_iterator iter = conf_clstr.begin();
                       iter != conf_clstr.end(); ++iter){
                               oarc << *iter;
                }

             
        }
     void load(graphlab::iarchive& iarc) {
               
               std::string v,vk;
               size_t num = 0, f=0;
               
             
               iarc >> commID;

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
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
                  //if(vk != "interest")
                  //std::cout <<"key: "<<vk<<", val: "<<v<< "\n";
               }

               iarc >> num;
              
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk ;
                  att_group.push_back(vk);
               }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> f;
                  freq.push_back(f);
               }
               double d;
               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> d;
                  CAG_supp.push_back(d);
               }
               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> d;
                  CAG_conf.push_back(d);
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
                  iarc >> f;
                  size_clstr.push_back(f);
               }
               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> d;
                  supp_clstr.push_back(d);
               }
               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> d;
                  conf_clstr.push_back(d);
               }

        }
 };


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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//distributed graph definition
typedef graphlab::distributed_graph<vdata, edata> graph_type;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#################################### Send data Items part #####################################
//

struct data_gatherer_directN {

        std::multimap<std::string, std::string> att_values;
        std::vector< std::multimap<std::string, std::string> > frnds_profiles;


       data_gatherer_directN () : att_values() {}
        void save(graphlab::oarchive& oarc) const {
           size_t num = frnds_profiles.size();
                oarc << num;
            for(size_t c=0; c<num; c++){
                oarc <<frnds_profiles[c].size();
                for(std::multimap<std::string, std::string>::const_iterator iter = frnds_profiles[c].begin();
                       iter != frnds_profiles[c].end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                      std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
                }
            }
        }
        void load(graphlab::iarchive& iarc) {
             size_t num, count;
              frnds_profiles.clear();
             std::string vk, v;
             iarc >> count;
             for(int j=0; j<count;j++){
             att_values.clear();

             iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
               }
            frnds_profiles.push_back(att_values);
           }
        }
   
       data_gatherer_directN& operator+=(const data_gatherer_directN& other) {
              for(std::vector<std::multimap<std::string, std::string> >::const_iterator iter = other.frnds_profiles.begin();
                       iter != other.frnds_profiles.end(); ++iter){
                  frnds_profiles.push_back(*iter);
                   //if(iter->first!="interest")
                   //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
              }

              return *this;
         }

};

class collect_attributes:public graphlab::ivertex_program<graph_type,data_gatherer_directN,data_message>, public graphlab::IS_POD_TYPE {

public:
 
       std::multimap<std::string, std::string> att_values;
       std::vector< std::multimap<std::string, std::string> > frnds_profiles;

       std::vector<single_message> msg2frwd;
       std::vector<std::string> att_group;
       std::vector<size_t> freq;
       std::vector<double> CAG_supp;
       std::vector<double> CAG_conf;
       

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
             /*for(std::multimap<std::string, std::string>::const_iterator iter = vertex.data().att_values.begin();
                       iter != vertex.data().att_values.end(); ++iter) {
                       att_n_v.att_values.insert(std::make_pair(iter->first,iter->second));
                       //if( iter->first =="education")
                       //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
              }*/
              
            if(edge.source().data().att_values.size() >0){
               std::multimap<std::string, std::string> att_v;               
               for(std::multimap<std::string, std::string>::iterator iter = edge.source().data().att_values.begin(); 
                   iter!=edge.source().data().att_values.end(); iter++){
                   std::string vk = iter->first;
                   std::string v = iter->second;
                   att_v.insert(std::make_pair(vk, v));
                   //if( vk !="interest")
                   //std::cout <<"key: "<<vk<<", val: "<<v<< "\n";
              }
              att_n_v.frnds_profiles.push_back(att_v);
           }
          return att_n_v;
        }

        void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
        }


   //====================================================================================

   void init(icontext_type& context,
                const vertex_type& vertex, const data_message &msg) {
        //if(vertex.data().commID == vertex.id())  
        //context.cout() << "iteration: " << context.iteration() << " , number of active vertices: " << context.num_vertices()  <<"\n";
        if(context.iteration() >= 1 ){
          for(size_t i=0; i<msg.msgs.size() ; i++){
             single_message reply =  msg.msgs[i];
             //check if I'm central node of message CommID
             if(vertex.data().commID == vertex.id()){
                for(size_t x=0; x< reply.att_group.size(); x++)
                   att_group.push_back(reply.att_group[x]);
                for(size_t x=0; x< reply.freq.size(); x++)
                   freq.push_back(reply.freq[x]);
                for(size_t x=0; x< reply.supp.size(); x++)
                   CAG_supp.push_back(reply.supp[x]);
                for(size_t x=0; x< reply.conf.size(); x++)
                   CAG_conf.push_back(reply.conf[x]);

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

              for(size_t x =0; x<vertex.data().freq_clstr.size(); x++){
                  //vertex.data().supp_clstr[x] = double(vertex.data().freq_clstr[x])/double(vertex.data().supp_clstr[x]);
                  //vertex.data().supp_clstr[x] /=100;
                  //vertex.data().conf_clstr[x] = vertex.data().supp_clstr[x]/double(vertex.data().conf_clstr[x]);
                  //vertex.data().conf_clstr[x]/=100;
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
             //else
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
    std::vector<graphlab::vertex_id_type> frnds, c_frnds;
    std::string attribute, token, token_sub, token_list;
    bool  hasFrnds = false, hasAttrib=false;
    size_t cid;
    std::vector<std::string> att_names;
    std::multimap<std::string, std::string>  att_values,c_att_values;

    att_names.push_back("school");
    att_names.push_back("major");
    att_names.push_back("employer");
    att_names.push_back("places_lived");

/////////////////////////////////////////////////Load data    
    strm >> vid;
    strm >> cid;
    
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
  //if(cid == 61935 || vid == 61935 || vid==4851 || vid==438 || vid==698){
   graph.add_vertex(vid, vdata(frnds,att_names,att_values,0));
   graph.add_vertex(cid, vdata(c_frnds,att_names, c_att_values, 0));
    //if(isCent)
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
        std::vector<std::string> att_vec;
        std::vector<size_t> att_freq;
        if(v.data().att_group_clstr.size()>0){    
          strm << "============================================== \n";
          strm <<  v.id() << " , " << v.data().friends.size() << "\n";
          
          for(size_t i=0; i< v.data().att_group_clstr.size(); i++){
             double supp = double(v.data().freq_clstr[i])/double(v.data().supp_clstr[i]);
             double conf = supp/double(v.data().conf_clstr[i]);
             strm << v.data().att_group_clstr[i] << " , " << v.data().freq_clstr[i] << " , " << supp << " , " << conf << "\n";
          }
       }
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

   
   size_t MAX_ITERATIONS=9;
   clopts.get_engine_args().set_option("max_iterations", MAX_ITERATIONS);
   graphlab::omni_engine<collect_attributes> engine2(dc, graph, "sync", clopts);
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



