// Copyright 2018 The MathWorks, Inc.

#include <zmq.hpp>
#include <iostream>
#include <vector>
#include <utility>
#include <memory>

#include "statcal_util.hpp"
#include "acrobotclient.hpp"

namespace {

#define REQUEST_TIMEOUT 2500 //  msecs (> 1000)
#define REQUEST_RETRIES  3   //  Number of tries before we abandon

const char *DELIMITERS = " ,"; // <space> or ","

// class ZmqMgr for managing socket connection with the server
class ZmqMgr {
  public:
    ZmqMgr(const std::string &addr) : context(1), socket_addr(addr)
    {
    }

    ~ZmqMgr() {}

    void sendRequest(const std::string & request_str);

    void retrieveReply(zmq::message_t & reply, int retries_left = REQUEST_RETRIES);
    
  private:
    std::string socket_addr;
    zmq::context_t context;
    std::unique_ptr<zmq::socket_t> socket_ptr;

    std::unique_ptr<zmq::socket_t> createSocket();
};

// Helper function to send a request with input arguments to the server
void sendRequest_helper(void *zm, const double input_data[], unsigned int size, const unsigned int iter)
{
    std::vector<double> Vec;
    for (int i = 0; i < size; i++){
        Vec.push_back(input_data[i]);
    }
    
    Vec.push_back(static_cast<double>(iter));

    std::string request_str;
    encode_double_data(Vec, request_str);

    reinterpret_cast<ZmqMgr *>(zm)->sendRequest(request_str);
}

// Helper function to retrieve reply from the server and parse its results
std::pair<double, double> retrieveReply_helper(void *zm)
{
    zmq::message_t reply;

    reinterpret_cast<ZmqMgr *>(zm)->retrieveReply(reply);

    char *reply_str = static_cast<char*>(reply.data());
    
    std::vector<double> data;
    int len = decode_header(reply_str);
    decode_double_data(len, reply_str, data);

    return std::make_pair(data[0], data[1]);
}

// ZmqMgr class method sendRequest
void ZmqMgr::sendRequest(const std::string & request_str)
{
    if (!socket_ptr) {
        socket_ptr = createSocket();
    }
    zmq::message_t request;
    request.rebuild(request_str.size());
    memset(request.data(), 0, request_str.size());    
    memcpy(request.data (), request_str.c_str(), request_str.size());
    
    // std::cout << "Sending " << request_str << std::endl;
    socket_ptr->send(request);
}

// ZmqMgr class method retrieveReply
void ZmqMgr::retrieveReply(zmq::message_t & reply, int retries_left)
{
    assert(socket_ptr);
    
    while (retries_left) {
        //  Poll socket for a reply, with timeout
        zmq::pollitem_t items[] = { {*socket_ptr, 0, ZMQ_POLLIN, 0 } };
        zmq::poll (&items[0], 1, REQUEST_TIMEOUT);
        
        //  If we got a reply, process it
        if (items[0].revents & ZMQ_POLLIN) {
            
            socket_ptr->recv(&reply);
            char *reply_str = static_cast<char*>(reply.data());
            reply_str[reply.size()] = '\0';
            
            //std::cout << "Received: " << reply_str << std::endl;
            break;
        } else if (--retries_left == 0) {
            throw std::runtime_error("Server connection timed out");
        } else {
            std::cout << "No response from server, retrying … " << std::endl;
        }
    }
}

// ZmqMgr class method createSocket
std::unique_ptr<zmq::socket_t> ZmqMgr::createSocket()
{
    std::unique_ptr<zmq::socket_t> s_ptr(new zmq::socket_t(context, ZMQ_REQ));

    s_ptr->connect(socket_addr.c_str());
    int linger = 0;
    s_ptr->setsockopt (ZMQ_LINGER, &linger, sizeof (linger));
    // std::cout << "Connecting to stats calculator server" << std::endl;
        
    return s_ptr;
}

} // anonymous namespace

// Wrapper functions
void *setupruntimeresources_wrapper(const std::string & connStr)
{
    return reinterpret_cast<void *>(new ZmqMgr(connStr));
}

void start_wrapper(unsigned int *iter_ptr)
{
    *iter_ptr = 0;
}

void outputs_wrapper(void *zm, unsigned int *iter_ptr, double *y1_ptr, double *y2_ptr)
{
    if (*iter_ptr > 0) {
        double y1, y2;
        std::tie( y1, y2) = retrieveReply_helper(zm);
        
        *y1_ptr = y1;
        *y2_ptr = y2;

    } else {
         *y1_ptr = 0.0;
         *y2_ptr = 0.0;
    }
}

void update_wrapper(void *zm, unsigned int *iter_ptr, const double input_data[], const unsigned int size )
{
    (*iter_ptr)++;
    sendRequest_helper(zm, input_data, size, *iter_ptr);
}



void terminate_wrapper(void *zm, unsigned int *iter_ptr)
{
    // Retrieve the reply for the last request from mdlUpdate
    if (*iter_ptr > 0) {
        double y1, y2;
        std::tie(y1, y2) = retrieveReply_helper(zm);
        std::cout << "Last results:  " << y1 << " " << y2 << std::endl;
    }
}

void cleanupruntimeresouces_wrapper(void *zm)
{
    delete reinterpret_cast<ZmqMgr *>(zm);
}

