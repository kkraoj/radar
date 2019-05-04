#ifndef HELPERS_HH
#define HELPERS_HH

#include <system_error>
#include <vector>
#include <complex>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include <boost/align/aligned_allocator.hpp>

class unix_error : public std::system_error {
    std::string what_;

public:
    unix_error( const std::string & context )
        : system_error( errno, std::system_category() ),
          what_( context + ": " + std::system_error::what() )
    {}

    const char *what() const noexcept override { return what_.c_str(); }
};

inline int Check( const char * context, const int return_value )
{
  if ( return_value >= 0 ) { return return_value; }
  throw unix_error( context );
}

inline int Check( const std::string & context, const int return_value )
{
    return Check( context.c_str(), return_value );
}

class FileDescriptor
{
    int fd_;

public:
    FileDescriptor( const int fd ) : fd_( fd ) {}

    ~FileDescriptor() { Check( "close", ::close( fd_ ) ); }

    int fd_num() const { return fd_; }

    uint64_t size() const
    {
        struct stat file_info;
        Check( "fstat", fstat( fd_, &file_info ) );
        return file_info.st_size;
    }

    FileDescriptor( const FileDescriptor & other ) = delete;
    FileDescriptor & operator=( const FileDescriptor & other ) = delete;
};

class MMAP
{
    FileDescriptor fd_;
    uint64_t size_;
    void * data_;

public:
    MMAP( const std::string & filename )
        : fd_( Check( filename, open( filename.c_str(), O_RDONLY ) ) ),
          size_( fd_.size() ),
          data_( mmap( nullptr, size_, PROT_READ, MAP_SHARED, fd_.fd_num(), 0 ) )
    {
        if ( data_ == MAP_FAILED ) { throw unix_error( "mmap" ); }
    }

    ~MMAP() { Check( "munmap", munmap( data_, size_ ) ); }

    uint8_t operator[] ( const uint64_t index ) const
    {
        return *( static_cast<uint8_t *>( data_ ) + index );        
    }

    uint8_t at( const uint64_t index ) const {
	if ( index >= size_ ) { throw std::out_of_range( "MMAP::at(): " + std::to_string( index ) + " >= " + std::to_string( size_ ) ); }
	return operator[]( index );
    }
    
    uint64_t size() const { return size_; }

    MMAP( const MMAP & other ) = delete;
    MMAP & operator=( const MMAP & other ) = delete;
};

using Signal = std::vector<std::complex<float>,
			   boost::alignment::aligned_allocator<std::complex<float>, 64>>;

class DAT
{
    MMAP file_;

public:
    DAT( const std::string & filename )
        : file_( filename )
    {
        if ( file_.size() % 4 ) {
            throw std::runtime_error( "File size needs to be a multiple of 4 bytes." );
        }
    }

    uint64_t IQ_sample_count() const { return file_.size() / 4; }

    int16_t I( const uint64_t index ) const
    {
        return file_[ 4 * index ] | (file_[ 4 * index + 1 ] << 8);
    }

    int16_t Q( const uint64_t index ) const
    {
        return file_[ 4 * index + 2 ] | (file_[ 4 * index + 3 ] << 8);
    }

    void read( Signal & signal, const uint64_t offset ) const
    {
	const uint64_t breakpoint = std::min( signal.size(), IQ_sample_count() );
	
	for ( unsigned int index = offset; index < breakpoint; index++ ) {
	    signal[ index ] = { float( I( index ) ), float( Q( index ) ) };
	}

	/* zero fill to end */
	for ( unsigned int index = breakpoint; index < signal.size(); index++ ) {
	    signal[ index ] = 0;
	}
    }
};

#endif /* HELPERS_HH */
