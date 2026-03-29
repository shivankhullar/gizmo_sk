/* This is a generic code block designed for simple neighbor loops, so that they don't have to
be copy-pasted and can be generically optimized in a single place */
{
    int j, k, ndone=0, ndone_flag=0, recvTask, place, save_NextParticle; long long n_exported = 0; double tstart, tend, tstart_loop; /* define some variables used only below */
    NextParticle = 0;    /* begin the main loop; start with this index into ActiveParticleList */
    tstart_loop = my_second();
    memset(ProcessedFlag, 0, All.MaxPart * sizeof(unsigned char));
    BufferCollisionFlag = 0; /* set to zero before operations begin */
    do /* primary point-element loop */
    {
        BufferFullFlag = 0; Nexport = 0; save_NextParticle = NextParticle; tstart = my_second();
        for(j = 0; j < NTask; j++) {Send_count[j] = 0; Exportflag[j] = -1;} /* do local particles and prepare export list */
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            PRIMARY_SUBFUN_NAME(&mainthreadid, loop_iteration);    /* do local particles and prepare export list */
        }
        tend = my_second(); timecomp += timediff(tstart, tend);
        if(BufferFullFlag) /* we've filled the buffer or reached the end of the list, prepare for communications */
        {
            int last_nextparticle = NextParticle;
            int processed_particles = 0;
            int first_unprocessedparticle = -1;
            NextParticle = save_NextParticle; /* figure out where we are */
            while(NextParticle < (int)ActiveParticleList.size())
            {
                if(NextParticle == last_nextparticle) {break;}
                int pindex = ActiveParticleList[NextParticle];
#ifndef _OPENMP
                if(ProcessedFlag[pindex] != 1) {break;}
#else
                if(ProcessedFlag[pindex] == 0 && first_unprocessedparticle < 0) {first_unprocessedparticle = NextParticle;}
                if(ProcessedFlag[pindex] == 1)
#endif
                {
                    processed_particles++;
                    ProcessedFlag[pindex] = 2;
                }
                NextParticle++;
            }
#ifdef _OPENMP
            if(first_unprocessedparticle >= 0) {NextParticle = first_unprocessedparticle;} /* reset the neighbor list properly for the next group since we can get 'jumps' with openmp active */
            if(processed_particles == 0 && NextParticle == save_NextParticle && NextParticle < (int)ActiveParticleList.size()) {
                BufferCollisionFlag++; if(BufferCollisionFlag < 2) {continue;}} /* we overflowed without processing a single particle, but this could be because of a collision, try once with the serialized approach, but if it fails then, we're truly stuck */
            else if(processed_particles && BufferCollisionFlag) {BufferCollisionFlag = 0;} /* we had a problem in a previous iteration but things worked, reset to normal operations */
#endif
            if(processed_particles <= 0 && NextParticle == save_NextParticle) // this is still sometimes being triggered with OPENMP, but not without, when it shouldn't. some OPENMP error still needs to be debugged
            {
                int pindex_dbg = (NextParticle < (int)ActiveParticleList.size()) ? ActiveParticleList[NextParticle] : -1;
                PRINT_WARNING("NextParticle == save_NextParticle condition (the buffer appears too small to hold a single particle): NextParticle=%d save_NextParticle=%d last_nextparticle=%d ProcessedFlag=%d NumPart=%d N_gas=%d NTaskTimesNumPart=%llu maxThreads=%d All.BunchSize=%ld All.BufferSize=%llu Nexport=%ld ndone=%d ndone_flag=%d NTask=%d",NextParticle,save_NextParticle,last_nextparticle,(pindex_dbg>=0 ? (int)ProcessedFlag[pindex_dbg] : -1),NumPart,N_gas,(unsigned long long)NTaskTimesNumPart,maxThreads,All.BunchSize,(unsigned long long)All.BufferSize,Nexport,ndone,ndone_flag,NTask);
                if(pindex_dbg >= 0) {PRINT_WARNING("This is a live particle: pindex=%d ID=%llu Mass=%g Type=%d",pindex_dbg,(unsigned long long)P[pindex_dbg].ID,P[pindex_dbg].Mass,P[pindex_dbg].Type);}
                endrun(113312);
            } /* in this case, the buffer is too small to process even a single particle */

            int new_export = 0; /* actually calculate exports [so we can tell other tasks] */
            for(j = 0, k = 0; j < Nexport; j++)
            {
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1) {k = j + 1;}
                    for(; k < Nexport; k++)
                        if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                        {
                            int old_index = DataIndexTable[j].Index;
                            DataIndexTable[j] = DataIndexTable[k]; DataNodeList[j] = DataNodeList[k]; DataIndexTable[j].IndexGet = j; new_export++;
                            DataIndexTable[k].Index = old_index; k++;
                            break;
                        }
                }
                else {new_export++;}
            }
            Nexport = new_export; /* counting exports... */
        }
        n_exported += Nexport;
        for(j = 0; j < NTask; j++) {Send_count[j] = 0;}
        for(j = 0; j < Nexport; j++) {Send_count[DataIndexTable[j].Task]++;}
        mysort_dataindex(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare); /* construct export count tables */
        tstart = my_second();
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD); /* broadcast import/export counts */
        tend = my_second(); timewait += timediff(tstart, tend);

        for(j = 0, Send_offset[0] = 0; j < NTask; j++) {if(j > 0) {Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];}} /* calculate export table offsets */
        DATAIN_NAME = (struct INPUT_STRUCT_NAME *) mymalloc("DATAIN_NAME", Nexport * sizeof(struct INPUT_STRUCT_NAME));
        DATAOUT_NAME = (struct OUTPUT_STRUCT_NAME *) mymalloc("DATAOUT_NAME", Nexport * sizeof(struct OUTPUT_STRUCT_NAME));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(j = 0; j < Nexport; j++) /* prepare particle data for export [fill in the structures to be passed] */
        {
            int place_local = DataIndexTable[j].Index;
            INPUTFUNCTION_NAME(&DATAIN_NAME[j], place_local, loop_iteration);
            memcpy(DATAIN_NAME[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }

        /* ok now we have to figure out if there is enough memory to handle all the tasks sending us their data, and if not, break it into sub-chunks */
        int N_chunks_for_import, ngrp_initial, ngrp;
        int max_ngrp = (1 << PTask);
        MPI_Request *mpi_reqs_import = (MPI_Request *) malloc(max_ngrp * sizeof(MPI_Request));
        MPI_Request *mpi_reqs_export = (MPI_Request *) malloc(max_ngrp * sizeof(MPI_Request));
        MPI_Request *mpi_reqs_result = (MPI_Request *) malloc(2 * max_ngrp * sizeof(MPI_Request));
        for(ngrp_initial = 1; ngrp_initial < (1 << PTask); ngrp_initial += N_chunks_for_import) /* sub-chunking loop opener */
        {
            int flagall;
            N_chunks_for_import = (1 << PTask) - ngrp_initial;
            do {
                int flag = 0; Nimport = 0;
                for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++)
                {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask) {if(Recv_count[recvTask] > 0) {Nimport += Recv_count[recvTask];}}
                }
                size_t space_needed = Nimport * sizeof(struct INPUT_STRUCT_NAME) + Nimport * sizeof(struct OUTPUT_STRUCT_NAME) + 16384; /* extra bitflag is a padding, to avoid overflows */
                if(space_needed > FreeBytes) {flag = 1;}

                MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                if(flagall) {N_chunks_for_import /= 2;} else {break;}
            } while(N_chunks_for_import > 0);
            if(N_chunks_for_import == 0) {printf("Memory is insufficient for even one import-chunk: N_chunks_for_import=%d  ngrp_initial=%d  Nimport=%ld  FreeBytes=%lld , but we need to allocate=%lld \n",N_chunks_for_import, ngrp_initial, Nimport, (long long)FreeBytes,(long long)(Nimport * sizeof(struct INPUT_STRUCT_NAME) + Nimport * sizeof(struct OUTPUT_STRUCT_NAME) + 16384)); endrun(9977);}
            if(flagall) {if(ThisTask==0) PRINT_WARNING("Splitting import operation into sub-chunks as we are hitting memory limits (check this isn't imposing large communication cost)");}

            /* now allocate the import and results buffers */
            DATAGET_NAME = (struct INPUT_STRUCT_NAME *) mymalloc("DATAGET_NAME", Nimport * sizeof(struct INPUT_STRUCT_NAME));
            DATARESULT_NAME = (struct OUTPUT_STRUCT_NAME *) mymalloc("DATARESULT_NAME", Nimport * sizeof(struct OUTPUT_STRUCT_NAME));

            /* Pipelined communication-computation overlap: split the ngrp range into two halves
               balanced by import count. Post all Irecv/Isend upfront, then process each half:
               wait for half1 imports → compute half1 (while half2 imports arrive in background)
               → post half1 results → wait for half2 → compute half2 (while half1 results exchange) */
            long Nimport_total = Nimport, Nimport_half1 = 0;
            int half_boundary = ngrp_initial + N_chunks_for_import;
            {
                long half_target = Nimport_total / 2, running = 0;
                for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++) {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask) {running += Recv_count[recvTask];}
                    if(running >= half_target) {half_boundary = ngrp + 1; Nimport_half1 = running; break;}
                }
                if(half_boundary == ngrp_initial + N_chunks_for_import) {Nimport_half1 = Nimport_total;}
            }

            /* Post ALL import Irecv's (both halves) and ALL export Isend's */
            tstart = my_second();
            int n_irecv = 0, n_irecv_half1 = 0, n_isend_export = 0;
            Nimport = 0;
            for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++)
            {
                recvTask = ThisTask ^ ngrp;
                if(recvTask < NTask)
                {
                    if(Recv_count[recvTask] > 0) {
                        MPI_Irecv(&DATAGET_NAME[Nimport], Recv_count[recvTask] * sizeof(struct INPUT_STRUCT_NAME),
                                  MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_A,
                                  MPI_COMM_WORLD, &mpi_reqs_import[n_irecv++]);
                    }
                    if(Send_count[recvTask] > 0) {
                        MPI_Isend(&DATAIN_NAME[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct INPUT_STRUCT_NAME),
                                  MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_A,
                                  MPI_COMM_WORLD, &mpi_reqs_export[n_isend_export++]);
                    }
                    Nimport += Recv_count[recvTask];
                }
                if(ngrp + 1 == half_boundary) {n_irecv_half1 = n_irecv;}
            }

            /* --- Half 1: wait for imports, compute, post results --- */
            MPI_Waitall(n_irecv_half1, mpi_reqs_import, MPI_STATUSES_IGNORE);
            tend = my_second(); timecomm += timediff(tstart, tend);

            tstart = my_second(); NextJ = 0; Nimport = Nimport_half1;
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
#ifdef _OPENMP
                int mainthreadid = omp_get_thread_num();
#else
                int mainthreadid = 0;
#endif
                SECONDARY_SUBFUN_NAME(&mainthreadid, loop_iteration);
            }
            tend = my_second(); timecomp += timediff(tstart, tend);

            int n_result_reqs = 0;
            {
                long nimport_offset = 0;
                for(ngrp = ngrp_initial; ngrp < half_boundary; ngrp++) {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask) {
                        if(Recv_count[recvTask] > 0) {
                            MPI_Isend(&DATARESULT_NAME[nimport_offset], Recv_count[recvTask] * sizeof(struct OUTPUT_STRUCT_NAME),
                                      MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_B,
                                      MPI_COMM_WORLD, &mpi_reqs_result[n_result_reqs++]);
                        }
                        if(Send_count[recvTask] > 0) {
                            MPI_Irecv(&DATAOUT_NAME[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct OUTPUT_STRUCT_NAME),
                                      MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_B,
                                      MPI_COMM_WORLD, &mpi_reqs_result[n_result_reqs++]);
                        }
                        nimport_offset += Recv_count[recvTask];
                    }
                }
            }

            /* --- Half 2: wait for imports, compute, post results --- */
            tstart = my_second();
            MPI_Waitall(n_irecv - n_irecv_half1, &mpi_reqs_import[n_irecv_half1], MPI_STATUSES_IGNORE);
            tend = my_second(); timecomm += timediff(tstart, tend);

            tstart = my_second(); NextJ = Nimport_half1; Nimport = Nimport_total;
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
#ifdef _OPENMP
                int mainthreadid = omp_get_thread_num();
#else
                int mainthreadid = 0;
#endif
                SECONDARY_SUBFUN_NAME(&mainthreadid, loop_iteration);
            }
            tend = my_second(); timecomp += timediff(tstart, tend);

            {
                long nimport_offset = Nimport_half1;
                for(ngrp = half_boundary; ngrp < ngrp_initial + N_chunks_for_import; ngrp++) {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask) {
                        if(Recv_count[recvTask] > 0) {
                            MPI_Isend(&DATARESULT_NAME[nimport_offset], Recv_count[recvTask] * sizeof(struct OUTPUT_STRUCT_NAME),
                                      MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_B,
                                      MPI_COMM_WORLD, &mpi_reqs_result[n_result_reqs++]);
                        }
                        if(Send_count[recvTask] > 0) {
                            MPI_Irecv(&DATAOUT_NAME[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct OUTPUT_STRUCT_NAME),
                                      MPI_BYTE, recvTask, TAG_MPI_GENERIC_COM_BUFFER_B,
                                      MPI_COMM_WORLD, &mpi_reqs_result[n_result_reqs++]);
                        }
                        nimport_offset += Recv_count[recvTask];
                    }
                }
            }

            /* Wait for all result exchanges and export sends to complete */
            tstart = my_second();
            MPI_Waitall(n_result_reqs, mpi_reqs_result, MPI_STATUSES_IGNORE);
            MPI_Waitall(n_isend_export, mpi_reqs_export, MPI_STATUSES_IGNORE);
            tend = my_second(); timecomm += timediff(tstart, tend);
            myfree(DATARESULT_NAME); myfree(DATAGET_NAME);

        } /* close the sub-chunking loop: for(ngrp_initial = 1; ngrp_initial < (1 << PTask); ngrp_initial += N_chunks_for_import) */
        free(mpi_reqs_result); free(mpi_reqs_export); free(mpi_reqs_import);

        /* we have all our results back from the elements we exported: add the result to the local elements */
        tstart = my_second();
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            OUTPUTFUNCTION_NAME(&DATAOUT_NAME[j], place, 1, loop_iteration);
        }
        tend = my_second(); timecomp += timediff(tstart, tend);
        myfree(DATAOUT_NAME); myfree(DATAIN_NAME); /* free the structures used to prepare our initial export data, we're done here! */

        if(NextParticle >= (int)ActiveParticleList.size()) {ndone_flag = 1;} else {ndone_flag = 0;} /* figure out if we are done with the particular active set here */
        tstart = my_second();
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /* call an allreduce to figure out if all tasks are also done here, otherwise we need to iterate */
        tend = my_second(); timewait += timediff(tstart, tend);
    }
    while(ndone < NTask);
    timeall += timediff(tstart_loop, my_second());

} /* closes clause, so variables don't 'leak' */
