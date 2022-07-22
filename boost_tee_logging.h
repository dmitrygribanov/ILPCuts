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
        string _log_file_name;
        ofstream _file;
        tee_device _tee_dev;
        stream _stream;

    public:
        logger(string log_file_name = "log.txt", ostream &a_ostr = cout)
            : _log_file_name(log_file_name), _file(_log_file_name), _tee_dev(a_ostr, _file), _stream(_tee_dev)
        {
        }

        string get_log_file_name() const { return _log_file_name; }

        void reset_output_file(string log_file_name)
        {
            _file.close();
            _file.open(log_file_name);
        }

        stream &out() { return _stream; }

        void print_command_line_params(int argc, char **argv)
        {
            out() << "ILPCuts command line params string:\n";
            for (int i = 0; i < argc; ++i)
            {
                out() << argv[i] << " ";
            }
            out() << "\n"
                  << endl;
        }
    };

} // namespace ILPCuts
