#pragma once

#include <ostream>
#include <fstream>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>

namespace ILPCuts
{
    class logger
    {
    public:
        using tee_device = boost::iostreams::tee_device<ostream, ofstream>;
        using stream = boost::iostreams::stream<tee_device>;

    private:
        ofstream _file;
        tee_device _tee_dev;
        stream _stream;

    public:
        logger(string log_file_name = "log.txt", ostream &a_ostr = cout)
            : _file(log_file_name), _tee_dev(a_ostr, _file), _stream(_tee_dev)
        {
        }

        stream &out() { return _stream; }
    };

} // namespace ILPCuts
