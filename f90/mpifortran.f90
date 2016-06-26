program mpi
	include 'mpif.h'
	
TYPE LIMITS
	integer start  
	integer stop
END TYPE

	integer, parameter :: MASTER = 0

integer i, j, k, N, procs, start, finish, rank, num_tasks, num_workers
integer status(MPI_STATUS_SIZE), threads, numToDo(*), len, dim, numRows
type(LIMITS) variableIndex(*)
double precision tmp(*), start_time, stop_time, trace, A(*), B(*), C(*)

pointer(ptr_A, A)
pointer(ptr_B, B)
pointer(ptr_C, C)
pointer(ptr_numToDo, numToDo)
pointer(ptr_variableIndex, variableIndex)
pointer(ptr_tmp, tmp)

	call MPI_INIT( ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, num_tasks, ierr ) 
	call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr ) 
	
	threads = num_tasks
	N = 5000
	len = N
	dim = N
	
ptr_A = malloc(N*N*8)
ptr_B = malloc(N*N*8)
ptr_C = malloc(N*N*8)
ptr_numToDo = malloc(threads*4)
ptr_variableIndex = malloc(threads*8)
	
	!Divide the code equally
	if (mod(len, threads) .eq. 0) then
		do i = 1, threads
			numToDo(i) = len / threads
		enddo
	else
		do i = 1, threads
			numToDo(i) = len / threads
		enddo
		do i = 1, mod(len, threads)
			numToDo(i) = numToDo(i) + 1
		enddo
	endif
	
	start = 0
	do i = 1, threads
		variableIndex(i)%start = start + 1
		variableIndex(i)%stop = start + numToDo(i)
		start = variableIndex(i)%stop
	enddo
	
	!Master function
	
	if (rank .eq. MASTER) then
		
		!Fill the matrices
		do i = 0, N - 1
			do j = 1, N
				A(i*N+j) = 1.0/sqrt(dble(N))
				B(i*N+j) = 1.0/sqrt(dble(N))
				C(i*N+j) = 0
			enddo
		enddo
		
		start_time = MPI_Wtime()
		
		!SEND DATA TO EVERY MACHINE
        !This code currently sends both matrices to all the machines.  The transfer
        !time cound be cut down significantly by only sending the rows of the A matrix
        !or the columns of the B matrix.
				do i = 1, threads - 1
			print *, "Sending A and B to process ", i
			call MPI_SEND(A, dim*dim, MPI_DOUBLE, i, MASTER, MPI_COMM_WORLD) !Send matrix A
			call MPI_SEND(B, dim*dim, MPI_DOUBLE, i, MASTER, MPI_COMM_WORLD) !Send matrix B
			print *, "Complete!"
		enddo
	
		!DO SECTION FOR WHICH MASTER IS RESPONSIBLE
		start = variableIndex(rank + 1)%start
		finish = variableIndex(rank + 1)%stop
		numRows = finish - start + 1
		
		ptr_tmp = malloc(dim*numRows*8)
		
		do i = 1, dim*numRows
			tmp(i) = 0.0
		enddo
		
		do i = start, finish
			do j = 1, dim
				do k = 1, dim
					tmp((i-start)*dim+j) = tmp((i-start)*dim+j) + A((i-1)*dim+k) * B((k-1)*dim+j)
				enddo
			enddo
		enddo
		
		do i = start, finish
			do j = 1, dim
				C((i-1)*dim+j) = tmp((i-start)*dim+j)
			enddo
		enddo
		
		call free(ptr_tmp)
		
		
		!RECEIVE DATA FROM EVERY MACHINE		
		do k = 1, threads - 1
			start = variableIndex(k + 1)%start
			finish = variableIndex(k + 1)%stop
			numRows = finish - start + 1
			
			ptr_tmp = malloc(dim*numRows*8)
			call MPI_RECV(tmp, dim*numToDo(k), MPI_DOUBLE, k, 0, MPI_COMM_WORLD, status)
			
			do i = start, finish
				do j = 1, dim
					C((i-1)*dim+j) = tmp((i-start)*dim+j)
				enddo
			enddo
			call free(ptr_tmp)	
		enddo
		
		stop_time = MPI_Wtime()
		trace = 0.0
		
		do i = 1, dim
			trace = trace + C((i-1)*dim+i)
		enddo
		
		print *, "Dimension: ", dim, " Trace: ", trace, " in ", stop_time-start_time, " seconds"	
	else
		!Worker Function ---------------------------------------------------
		
		call MPI_RECV(A, dim*dim, MPI_DOUBLE, MASTER, MASTER, MPI_COMM_WORLD, status)
		call MPI_RECV(B, dim*dim, MPI_DOUBLE, MASTER, MASTER, MPI_COMM_WORLD, status)
		
		start = variableIndex(rank + 1)%start
		finish = variableIndex(rank + 1)%stop
		numRows = finish - start + 1
		
		ptr_tmp = malloc(dim*numRows*8)
		do i = 1, dim*numRows
			tmp(i) = 0.0
		enddo
		
		do i = start, finish
			do j = 1, dim
				do k = 1, dim
					tmp((i-start)*dim+j) = tmp((i-start)*dim+j) + A((i-1)*dim+k) * B((k-1)*dim+j)
				enddo
			enddo
		enddo
		
		print *, "Worker ", rank, " sending to Master"
		call MPI_SEND(tmp, dim*numRows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD)
		print *, "Worker ", rank, " finished sending"
		
		call free(ptr_tmp)
	endif
	
	call MPI_FINALIZE()
end program
