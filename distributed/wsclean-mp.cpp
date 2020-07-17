#include "slave.h"

#include "taskmessage.h"

#include "../wsclean/commandline.h"
#include "../wsclean/logger.h"
#include "../wsclean/wsclean.h"

#include <exception>
#include <iostream>

#include <mpi.h>

int main(int argc, char *argv[])
{
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided != MPI_THREAD_MULTIPLE)
	{
		std::cout << "This MPI implementation does not support multiple threads.\n";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	int result = 0;
	WSClean wsclean;
	int world_size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	bool master = (rank == 0);
	
	char hostname[256];
	gethostname(hostname, sizeof(hostname));
	std::cout << "Node " << rank << ", PID " << getpid() << " on " << hostname << "\n";
	
	bool shortException = false;
	try {
		bool parseResult = false;
		shortException = !master;
		parseResult = CommandLine::Parse(wsclean, argc, argv, !master);
		shortException = false;
		if(parseResult)
		{
			WSCleanSettings& settings = wsclean.Settings();
			settings.useMPI = true;
			
			if(master) {
				CommandLine::Run(wsclean);
				TaskMessage message;
				message.type = TaskMessage::Finish;
				for(int i=1; i!=world_size; ++i)
				{
					MPI_Send(
						&message,
						sizeof(TaskMessage),
						MPI_BYTE,
						i,
						0,
						MPI_COMM_WORLD);
				}
			}
			else {
				Slave slave(settings);
				slave.Run();
			}
		}
		Logger::Error << "Process " << rank << " finished.\n";
	} catch(std::exception& e)
	{
		if(shortException)
			Logger::Error << "Process " << rank << " stopped because of an exception.\n";
		else {
			Logger::Error
				<< "+ + + + + + + + + + + + + + + + + + +\n"
				<< "+ An exception occured in process " << rank << ":\n"
				<< "+ >>> " << e.what() << "\n"
				<< "+ + + + + + + + + + + + + + + + + + +\n";
		}
		result = -1;
	}
	MPI_Finalize();
	return result;
}

