/*--------------------------------------------------------------------*\
 * file -- task_prepare.c
 * V-2.0_0
 *
 * function -- int task_prepare(List * tasks)
 * - preparing new tasks in the presence of idle CPUs;
 *
 * (WARNING!) NOT-ADAPTIVE ALGORITHM - only checks the number of
 * free CPUs of the selected mashine at a time without taking into
 * account the presence of processes in the state of waiting-switching
 * to execution. May be brakes.
 *
 * 1. get the number of free CPUs;
 * 2. rank mashines with free cpu by performance;
 *
 * ENTER to CYCLE up to data-sources file exists AND the counter
 * of prepared tasks LT the number of free CPUs...
 *
 * 3. choose one of the most productive machines with free cpu;
 * 4. choose the largest sources-data file in the struct `tasks_list'
 *    marked as 0 - waiting in `status' field;
 * 5. create a working directory (if not exists) for the selected
 *    machine/task-name in the shared area of local FS, named as
 *    `TASK_NAME'/.IP_ADDR_bind_name_of_blade_selected/ included:
 *    {work/,work/sources/,work/out}, further as a '.../ in the text';
 * 6. cooy ---> to: - program ---> .../work/ if not exists,
 *                  - sources-data file --> .../work/sources/;
 * 7. mark as 1 - `prepared' in `struct tasks_list';
 * 8. increase the counter of prepared tasks;
 * 9. if data-sources file not exists ?
 *
 * EXIT from Cycle;
 *
 * 10. Return counter of prepared tasks, the end.
 *
 * WHERE:
 *  tasks - pointer to the head of linked `List' which include elems.
 *          of TYPE `Item' described an infomation of each will
 *          prepared to launch tasks.
 *
 * RETURNS: 1...N - tasks (source data file) were prepared and may be
 *                  started;
 *          0 - no one task (source data file) were prepared - this
 *              may be if there are no free CPUs or all tasks are
 *              in progress;
 *         -1 - error.
\*--------------------------------------------------------------------*/
static int x, y, z;
static int get_maxsrc_cutfile();

//######################################################################
//----------------------------------------------------------------------
int
task_prepare(List * tasks)
{
  int ret_code = 0;
  int ret_func;
  // counters:
  int i, j, k, n;
  int count_tasks_prep = 0;
  int count_freecpu = 0;
  int count_tot_freecpu = 0;
  int max_strsrc = 0;                                                   // the returned value by the func. get_maxsrc_cutfile();
  // flags:
  int flag_remotefunc;                                                  // 0 - if remote proc. returns an none successful-status;
  int flag_chfreq;                                                      // if nodes of linked lists have replaced;
  int flag_freecpu;
  // handels:
  extern CLIENT * handle;                                               // handle of the established connection with a remote mashine (with remote procs.);
  // structures & pointers:
  extern struct blade_tasks tasks_list[MAX_TASKS];                      //
  extern struct calc_hosts super_hosts[MAX_BLADES];                     //
  Item item_temp;                                                       // create a new Item (not included into the List of nodes), will copy to list;
  Node * pscan_tasks;                                                   // pointer to the top of the list nodes;
  // conf params:
  extern char addr0[ADDRLN];                                            // at control-host interface InfiniBand Factory Over IP;
  extern char addr1[ADDRLN];                                            // at control-host interface Ethernet;


  char pathname_machinename[PATHLN];
  char pathname_taskname[PATHLN];
  char pathname_task[PATHLN*2];
  char save_taskname[PATHLN];
  char shcommand[PATHLN*2];
  char testfile[PATHLN];

  // get system parameters of the remote machines and save its in to the struct. named `superhosts';
  for (i = 0; i < MAX_BLADES; i++)
  {
    if (strlen(super_hosts[i].chassis_name) > 0 && super_hosts[i].online == 1)
    {
      if (strcmp(super_hosts[i].addr_eth, "127.0.0.1") == 0 || strcmp(super_hosts[i].addr_eth, addr0) == 0 || strcmp(super_hosts[i].addr_eth, addr1) == 0) // if local machine;
      {
        super_hosts[i].cpu_numb = local_getnumbcpu();
        super_hosts[i].cpu_freq = local_getfreqcpu();
        super_hosts[i].mem_total = local_gettotmem();
        super_hosts[i].mem_tasks = local_getmaxprocmem();
        super_hosts[i].procs_active = local_getnumbrun();
        super_hosts[i].manager_chez = (int)local_alreadyrun("__multi_CHLD");
        if (super_hosts[i].manager_chez > 0)
          super_hosts[i].procs_chez = local_getnumbchld((pid_t)super_hosts[i].manager_chez);
        continue;
      }
      // if a remote machine;
      //##### ESTABLISHING A CONNECTION TO CALL A REMOTE PROCEDURE #####
      handle = clnt_create(super_hosts[i].addr_eth, RSUPERBLPROG, RSUPERBLVERS, "tcp"); // client creation routine for program `prognum' and `version';
      if (handle == 0)
      {
        if (flag_testme == 1)                                                 // if self-test active;
          printf("could not contact remote program on mashine with interface ip address %s.\n", super_hosts[i].addr_eth);
        super_hosts[i].online = 0;                                    // switch-off in local structure;
        continue;
      }
      //######################################################################
      flag_remotefunc = 1;

      super_hosts[i].cpu_numb = hosts_getnumbcpu();
      if (super_hosts[i].cpu_numb <= 0) flag_remotefunc = 0;
      super_hosts[i].cpu_freq = hosts_getfreqcpu();
      if (super_hosts[i].cpu_freq <= 0) flag_remotefunc = 0;
      super_hosts[i].mem_total = hosts_gettotmem();
      if (super_hosts[i].mem_total <= 0) flag_remotefunc = 0;
      super_hosts[i].mem_tasks = hosts_getmaxprocmem();
      if (super_hosts[i].mem_tasks <= 0 || (super_hosts[i].mem_total - super_hosts[i].mem_tasks) <= 300) flag_remotefunc = 0;
/*
      super_hosts[i].procs_active = hosts_getnumbrun();
      if (super_hosts[i].procs_active >= super_hosts[i].cpu_numb) flag_remotefunc = 0;
*/
      // get number tasks have running;
      super_hosts[i].procs_active = 0;
      for (j = 0; j < MAX_TASKS; j++)
      {
        if (strlen(tasks_list[j].task_name) > 0)
        {
          for (k = 0; k < MAXSRCFILES; k++)
          {
            if (tasks_list[j].sources_items[k].numb_items > 0)
            {
              for (n = 0; n < MAXSRCFILESCUT; n++)
              {
                if (   tasks_list[j].sources_items[k].gen_source_items[n].status == 2
                    && tasks_list[j].sources_items[k].gen_source_items[n].calc_hosts_index == i
                    && hosts_ifpidexists(tasks_list[j].sources_items[k].gen_source_items[n].pid) == 1 )
                {
                  super_hosts[i].procs_active++;
                }
              }
            }
          }
        }
      }
/*
      super_hosts[i].manager_chez = hosts_alreadyrun("__multi_CHLD");
      if (super_hosts[i].manager_chez > 0)
        super_hosts[i].procs_chez = hosts_getnumbchld(super_hosts[i].manager_chez);
      if (super_hosts[i].manager_chez > 0 || super_hosts[i].procs_chez >= super_hosts[i].cpu_numb) flag_remotefunc = 0;
*/
      if (flag_remotefunc == 0)                                       // if have been an error...
        super_hosts[i].online = 0;                                    // ... mark remote mashane as off-line;
      clnt_destroy(handle);                                           // destroy client by the handle;
    }
    if (flag_testme == 1)
    {
      printf("task_prepare() #%d. {%s, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d, %d, %d, %d, %d}\n", i, super_hosts[i].chassis_name, \
                                                                                                      super_hosts[i].addr_eth, \
                                                                                                      super_hosts[i].addr_ib, \
                                                                                                      super_hosts[i].addr_elom, \
                                                                                                      super_hosts[i].acc, \
                                                                                                      super_hosts[i].pw, \
                                                                                                      super_hosts[i].bind_name, \
                                                                                                      super_hosts[i].cpu_numb, \
                                                                                                      super_hosts[i].cpu_freq, \
                                                                                                      super_hosts[i].mem_total, \
                                                                                                      super_hosts[i].mem_tasks, \
                                                                                                      super_hosts[i].procs_active, \
                                                                                                      super_hosts[i].manager_chez, \
                                                                                                      super_hosts[i].procs_chez, \
                                                                                                      super_hosts[i].online);
      fflush(stdout);
    }
  }

// 1. Get the number of free CPUs.

  // fill the structure with free_cpus;
  InitializeList(tasks);

  for (i=0; i<MAX_BLADES && strlen(super_hosts[i].chassis_name) > 0; i++)
  {
    if (super_hosts[i].online == 0)
      continue;
    count_freecpu = super_hosts[i].cpu_numb - super_hosts[i].procs_active - 2;
    if (count_freecpu > 0)
      count_tot_freecpu += count_freecpu;

    // add data in to the tmp struct;
    if (IsFullList(*tasks))                                              // check if a new node may be creted;
      errExit("IsFullList: can't do memory allocation a new item.");
    item_temp.index_host = i;
    item_temp.cpu_free = count_freecpu;
    item_temp.cpu_freq = super_hosts[i].cpu_freq;
    item_temp.sources_numb = 0;
    if (AppendItem(item_temp, tasks) == False)
      errExit("AppendItem: can't do memory allocation for a new item.\n");
  }

// 2. Rank mashines with free cpu by performance.
  i = 0;
  do                                                                    // A: find if current node Freq LT next and swap their;
  {
    flag_chfreq = 0;                                                    // 1 - if swap action have been;
    pscan_tasks = * tasks;                                              // set the top of the list;
    while (pscan_tasks != NULL)                                         // while not the end of listl
    {
      if (pscan_tasks->next != NULL)                                    // if current is not latest (have the next node);
      {
        if (pscan_tasks->item.cpu_freq < pscan_tasks->next->item.cpu_freq)  // if curren node-cpu-frequency LT next node-cpu-frequency...
        {
          if (SwapItemsUP(tasks, pscan_tasks->next))                    // ... swap the current and next nodes;
          {
            flag_chfreq = 1;                                            // ... if the replacement was successful...
            break;                                                      // ... interrupt cycle and go to the A:
          }
        }
      }
      pscan_tasks = pscan_tasks->next;                                  // go to the next node;
    }
  } while (flag_chfreq == 1 && i++ < 10000 );   // if rank successful or loop threat ... the end;

  if (flag_testme == 1)
  {
    print_list(*tasks);
    printf("\n\n LARGEST FILE:+++++++++++++++++++++++++++--%d--%d--%d--%d\n", x,y,z, get_maxsrc_cutfile());
  }

// ENTER to CYCLE until a waited data-sources file exists AND the counter of prepared tasks LT the number of free CPUs...
  count_tasks_prep = 0;
  pscan_tasks = * tasks;                                              // set the top of the list;
// 3. Choose one of the most productive machines with free cpu.
  while (pscan_tasks != NULL)                                         // while not the end of listl
  {
    flag_freecpu = 0;
    if (pscan_tasks->item.cpu_free > 0)   // if free-cpu exists... ;
    {
      flag_freecpu = 1;
// 4. Choose the largest sources-data file in the struct `tasks_list' marked as 0 - waiting in `status' field.
      if (get_maxsrc_cutfile() > 0)                                     // if waiting task exists... (sets x,y,z vector) ;
      {                                                                 // sets the x,y,z vector of a task as tasks_list[x].sources_items[y].gen_source_items[z] with max amount strings of waited tasks;
        //... CREATE RUN...
// 5. Create a working directory for the selected machine/task-name (if not exists) in the shared area of local FS, named as
//    `TASK_NAME'/.IP_ADDR_bind_name_of_blade_selected/ included: {work/,work/sources/,work/out}, further as a '.../ in the text'.

        // copy vector of task into the current node;
        pscan_tasks->item.sources[pscan_tasks->item.sources_numb].index_task    = x;
        pscan_tasks->item.sources[pscan_tasks->item.sources_numb].index_srcfile = y;
        pscan_tasks->item.sources[pscan_tasks->item.sources_numb].index_item    = z;

        // create task directory name;
        strncpy(pathname_task, tasks_list[x].task_name, PATHLN);
        for (n = 0; n < PATHLN && pathname_task[n] != '\0'; n++)        // to walk thru the string...
          if (isalpha(pathname_task[n]) == 0 && isdigit(pathname_task[n]) == 0)
            pathname_task[n] = '_';                                     // ... remove silo;

        snprintf( save_taskname, PATHLN, ".%s_%s/%s/work", \
                  super_hosts[pscan_tasks->item.index_host].addr_eth, \
                  super_hosts[pscan_tasks->item.index_host].bind_name, \
                  pathname_task );

        snprintf( pathname_machinename, PATHLN, "%s/.%s_%s", \
                  mount_of,
                  super_hosts[pscan_tasks->item.index_host].addr_eth, \
                  super_hosts[pscan_tasks->item.index_host].bind_name );

        snprintf( pathname_taskname, PATHLN, "%s/%s", \
                  pathname_machinename, \
                  pathname_task );

        snprintf( pathname_task, PATHLN, "%s/work", \
                  pathname_taskname );

        // create task directory by a machine;
        if (access(pathname_machinename, F_OK) != 0)
        {
          snprintf(shcommand, PATHLN*2, "mkdir -m 0777 %s", pathname_machinename);
          ret_func = system(shcommand);                                 // create ~ /nfs_share/work/.192.168.1.205_blade12one;
          if (ret_func < 0 || ret_func == 127)
          {
            if (flag_testme == 1)
              errMsg("WARNING: task_prepare() -- can't create a machine directory `%s'.", pathname_machinename);
            continue;
          }
        }
        // create task directory by a machine/task_name;
        if (access(pathname_taskname, F_OK) != 0)
        {
          snprintf(shcommand, PATHLN*2, "mkdir -m 0777 %s", pathname_taskname);
          ret_func = system(shcommand);                                   // create ~ /nfs_share/work/.192.168.1.205_blade12one/TASK_3;
          if (ret_func < 0 || ret_func == 127)
          {
            if (flag_testme == 1)
              errMsg("WARNING: task_prepare() -- can't create working directory `%s'.", pathname_taskname);
            continue;
          }
          // +/work as WR
          snprintf(shcommand, PATHLN*2, "mkdir -m 0777 %s", pathname_task);
          ret_func = system(shcommand);
          if (ret_func < 0 || ret_func == 127)
          {
            if (flag_testme == 1)
              errMsg("WARNING: task_prepare() -- can't create working directory `%s'.", pathname_task);
            continue;
          }
          // add +/sources RO
          snprintf(shcommand, PATHLN*2, "mkdir -m 0777 %s/sources", pathname_task);
          ret_func = system(shcommand);
          if (ret_func < 0 || ret_func == 127)
          {
            if (flag_testme == 1)
              errMsg("WARNING: task_prepare() -- can't create working directory `%s/sources'.", pathname_task);
            continue;
          }
          // add +/out as WR
          snprintf(shcommand, PATHLN*2, "mkdir -m 0777 %s/out", pathname_task);
          ret_func = system(shcommand);
          if (ret_func < 0 || ret_func == 127)
          {
            if (flag_testme == 1)
              errMsg("WARNING: task_prepare() -- can't create working directory `%s/out'.", pathname_task);
            continue;
          }
          // add +/configs as WR
          snprintf(shcommand, PATHLN*2, "mkdir -m 0777 %s/configs", pathname_task);
          ret_func = system(shcommand);
          if (ret_func < 0 || ret_func == 127)
          {
            if (flag_testme == 1)
              errMsg("WARNING: task_prepare() -- can't create working directory `%s/configs'.", pathname_task);
            continue;
          }
        }
//printf("\n+++++++++++++++++++++++++++++%s :: %d\n", shcommand, ret_func);
        // create source-dir and out-dir in the task directory;
// 6. Copy the task files...:
//      a) a program ---> .../work/ if not exists,
        snprintf(testfile, PATHLN, "%s/%s", pathname_task, strrchr(tasks_list[x].progr, '/'));
        if (access(testfile, F_OK) != 0)
        {
          snprintf(shcommand, PATHLN*2, "cp %s %s", tasks_list[x].progr, testfile);
          ret_func = system(shcommand);
          if (ret_func < 0 || ret_func == 127)
          {
            if (flag_testme == 1)
              errMsg("WARNING: task_prepare() -- can't copy file `%s' to `/%s'.", tasks_list[x].progr, testfile);
            continue;
          }
        }
//      b) sources-data file ---> .../work/sources/.
        snprintf(testfile, PATHLN, "%s/sources/", pathname_task);
        snprintf(shcommand, PATHLN*2, "cp %s %s/sources/", tasks_list[x].sources_items[y].gen_source_items[z].name, pathname_task);
        ret_func = system(shcommand);
        if (ret_func < 0 || ret_func == 127)
        {
          if (flag_testme == 1)
            errMsg("WARNING: task_prepare() -- can't copy src data file `%s' to `%s'.", tasks_list[x].sources_items[y].gen_source_items[z].name, testfile);
          continue;
        }
//      c) inversion configuration file ---> .../work/configs/.
        snprintf(testfile, PATHLN, "%s/configs/", pathname_task);
        snprintf(shcommand, PATHLN*2, "cp -R %s/ %s/configs/", "WORK/configs", pathname_task);
        ret_func = system(shcommand);
        if (ret_func < 0 || ret_func == 127)
        {
          if (flag_testme == 1)
            errMsg("WARNING: task_prepare() -- can't copy configs data file `%s' to `%s'.", tasks_list[x].sources_items[y].gen_source_items[z].name, testfile);
          continue;
        }

        // copy a task working directory name by remote mashine to `struct blade_tasks';
        strncpy(tasks_list[x].sources_items[y].gen_source_items[z].path_task, save_taskname, PATHLN);
// 7. Mark as 1 - `prepared' in `struct tasks_list'.
        tasks_list[x].sources_items[y].gen_source_items[z].status = 1;  // change the cut source data file status to `prepared';
        pscan_tasks->item.sources_numb++;                               // incr. the counter of sources data files of the current node;
        pscan_tasks->item.cpu_free--;                                   // decrement free-cpu common counter in the current list node;
// 8. Increase the counter of prepared tasks.
        count_tasks_prep++;                                             // increment tasks-prepared counter;

        if (flag_testme == 1)
        {
          // the data will rewrite to `struct blade_tasks tasks_list[MAX_TASKS]' after ruunig the task (proc. task_run());
          printf("--- %d +++++ %d \t %d \t %d \t [%d,%d,%d](%d) -- remote path:%s; source data file:%s; progr.name:%s;\n", \
                  count_tasks_prep, \
                  pscan_tasks->item.index_host, \
                  pscan_tasks->item.cpu_freq, \
                  pscan_tasks->item.cpu_free, \
                  x,y,z, pscan_tasks->item.sources_numb, \
                  tasks_list[x].sources_items[y].gen_source_items[z].path_task, \
                  strrchr(tasks_list[x].sources_items[y].gen_source_items[z].name, '/'), \
                  strrchr(tasks_list[x].progr, '/'));
        }
      }
// 9. If data-sources file not exists...
      else                                                              // WAITED TASK NOT FOUND ! ;
// EXIT from Cycle... (up to string ~350);
        break;                                                          // nothing more to do...

      if (pscan_tasks->next != NULL && pscan_tasks->next->item.cpu_freq < pscan_tasks->item.cpu_freq)  // if current is not latest (have the next node) and next node have LT freq... ;
      {
        pscan_tasks = * tasks;                                          // set pointer to the 1st node;
        continue;                                                       // run cycle from 1st node...;
      }
    }

    if (pscan_tasks->next == NULL && flag_freecpu > 0)
    {
      pscan_tasks = * tasks;                                          // set pointer to the 1st node;
      continue;                                                       // run cycle from 1-st node...;
    }
    else
      pscan_tasks = pscan_tasks->next;                                    // go to the next node;
  }

  // 10. Return counter of prepared tasks, the end.
  if (ret_code == 0)
    ret_code = count_tasks_prep;

  return ret_code;
}
/*####################################################################*/

/*--------------------------------------------------------------------*\
 * function - void get_maxsrc_cutfile()
 * V1.1_0
 * Stets point vector x,y,z to the run task (in struct blade_tasks)
 * with greatest number of strings.
 *
 * RETURN: 0...N -- the number of strings in the selected element
 *                  of the struct.
\*--------------------------------------------------------------------*/
static int
get_maxsrc_cutfile()
{
  extern int x,y,z;
  int i, j, k;
  int max_str = 0;

  for (i = 0; i < MAX_TASKS; i++)
  {
    if (strlen(tasks_list[i].task_name) > 0)
    {
      for (j = 0; j < MAXSRCFILES; j++)
      {
        if (tasks_list[i].sources_items[j].numb_items > 0)
        {
          for (k = 0; k < MAXSRCFILESCUT; k++)
          {
            if (tasks_list[i].sources_items[j].gen_source_items[k].numb_strings > 0 && tasks_list[i].sources_items[j].gen_source_items[k].status == 0)
            {
              if (tasks_list[i].sources_items[j].gen_source_items[k].numb_strings > max_str)
              {
                max_str = tasks_list[i].sources_items[j].gen_source_items[k].numb_strings;
                x = i;
                y = j;
                z = k;
//printf("+++++++++++++++++++++++++++--%d--%d--%d--%d\n", x,y,z, max_str);
              }
            }
          }
        }
      }
    }
  }

  return max_str;
}
/*####################################################################*/

/*####################################################################*/
/*##### BOTTOM OF THE FILE ###########################################*/

