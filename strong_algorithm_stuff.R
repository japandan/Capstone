#'A helper method that allows for multi-threading of calculating of strong peptide counts
#'This is an internal method used by process_strong and requires
#'the conn list of connection information to be passed into it by way of
#'cluster import.  This method will loop through all genus_id values in in_set
#'and populate a temporary table that stores the percentage of organisms in the
#'genus that reference a given peptide, for all peptides referenced by the genus.
#'
#'@param in_set a vector of genus_id values from the database
#'@param conn database connection values to create a MySQL connection
#'
#'@return the in_set vector that was passed into the method
#'
#'@importFrom DBI sqlInterpolate dbExecute dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'
#'@author Dustin Crockett
analyze_strong_counts <- function(in_set, conn) {
  
  #validate the conn list
  # verify that the conn is a list
  assertthat::assert_that(is.list(conn),
                          msg = "conn is not a list. Make sure conn is a list.")
  # verify that conn has the correct fields
  assertthat::assert_that("dbname" %in% names(conn),
                          msg = "dbname is missing from the list conn")
  assertthat::assert_that("host" %in% names(conn),
                          msg = "host is missing from the list conn")
  assertthat::assert_that("port" %in% names(conn),
                          msg = "port is missing from the list conn")
  assertthat::assert_that("user" %in% names(conn),
                          msg = "user is missing from the list conn")
  assertthat::assert_that("password" %in% names(conn),
                          msg = "password is missing from the list conn")
  
  con1 <- DBI::dbConnect(RMariaDB::MariaDB(),
                         dbname = conn$dbname,
                         host = conn$host,
                         port = conn$port,
                         user = conn$user,
                         password = conn$password)
  for(i in in_set) {
    sqlcmd <- DBI::sqlInterpolate(con1, "insert into candidate.strong_peptide_temp(peptide_id, genus_id, strong_percent)
                             select otp.peptide_id, s.genus_id, count(distinct o.id) / genus_count.org_count organism_percentage
                             from candidate.organisms o
                             join candidate.species s on s.id = o.species_id and s.genus_id = ?genus_id
                             join candidate.organisms_to_peptides otp on o.id = otp.organism_id
                             join (select count(distinct o.id) org_count, s.genus_id from candidate.organisms o join candidate.species s on s.id = o.species_id
                             where s.genus_id = ?genus_id
                             group by s.genus_id) genus_count
                             on s.genus_id = genus_count.genus_id
                             group by otp.peptide_id, s.genus_id, genus_count.org_count", genus_id = i)
    DBI::dbExecute(con1, sqlcmd)
  }
  DBI::dbDisconnect(con1)
  return(in_set)
}


#'A helper method that allows for multi-threading of calculating of strong peptide counts
#'This is an internal method used by process_strong and requires
#'the conn list of connection information to be passed into it by way of
#'cluster import.  This method will loop through all genus_id values in in_set
#'and populate a temporary table that stores the percentage of organisms in the
#'genus that reference a given peptide, for all peptides referenced by the genus.
#'
#'@param in_set a vector of genus_id values from the database
#'@param conn database connection values to create a MySQL connection
#'
#'@return the in_set vector that was passed into the method
#'
#'@importFrom DBI sqlInterpolate dbExecute dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'
#'@author Dustin Crockett
analyze_strong_counts <- function(in_set, conn) {
  
  #validate the conn list
  # verify that the conn is a list
  assertthat::assert_that(is.list(conn),
                          msg = "conn is not a list. Make sure conn is a list.")
  # verify that conn has the correct fields
  assertthat::assert_that("dbname" %in% names(conn),
                          msg = "dbname is missing from the list conn")
  assertthat::assert_that("host" %in% names(conn),
                          msg = "host is missing from the list conn")
  assertthat::assert_that("port" %in% names(conn),
                          msg = "port is missing from the list conn")
  assertthat::assert_that("user" %in% names(conn),
                          msg = "user is missing from the list conn")
  assertthat::assert_that("password" %in% names(conn),
                          msg = "password is missing from the list conn")
  
  con1 <- DBI::dbConnect(RMariaDB::MariaDB(),
                         dbname = conn$dbname,
                         host = conn$host,
                         port = conn$port,
                         user = conn$user,
                         password = conn$password)
  for(i in in_set) {
    sqlcmd <- DBI::sqlInterpolate(con1, "insert into candidate.strong_peptide_temp(peptide_id, genus_id, strong_percent)
                             select otp.peptide_id, s.genus_id, count(distinct o.id) / genus_count.org_count organism_percentage
                             from candidate.organisms o
                             join candidate.species s on s.id = o.species_id and s.genus_id = ?genus_id
                             join candidate.organisms_to_peptides otp on o.id = otp.organism_id
                             join (select count(distinct o.id) org_count, s.genus_id from candidate.organisms o join candidate.species s on s.id = o.species_id
                             where s.genus_id = ?genus_id
                             group by s.genus_id) genus_count
                             on s.genus_id = genus_count.genus_id
                             group by otp.peptide_id, s.genus_id, genus_count.org_count", genus_id = i)
    DBI::dbExecute(con1, sqlcmd)
  }
  DBI::dbDisconnect(con1)
  return(in_set)
}

#'Populate the database with strong peptides values
#'This method is multi threaded.  First creating a temporary table and then calling analyze_strong_counts
#'in parallel.  The number times analyze_strong_counts is called defined by the core_count variable,
#' the default value is 4.  Once analyze_strong_counts populates the
#'temporary table with the coverage of organisms by peptide and genus.  Next this temporary table
#'is agregated to create strong peptides records for all coverages that stand out from other coverages
#'by the tolerance provided in the tolerance variable, the
#'default value is .053.  Once that is done the temporary table is deleted.
#'
#'Throughout the population process information messages are printed out to show progress.
#'
#'`process_strong` calculates which peptides that are strong for a genus, for all peptides.
#'
#'
#'@param conn database connection values to create a MySQL connection
#'@param core_count This is to total number cores to be used during parallel processing
#'      it should be no greater than 1/2 of total cores
#'@param split_size This is size of list that should be passed to each execution of
#'      analyze_strong_count.
#'@param tolerance the percentage difference between organism coverages within a genus
#'      of the two highest covered genus, for the highest covered genus to be
#'      have a strong peptide
#'
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbExecute dbConnect dbDisconnect
#'@importFrom RMariaDB MariaDB
#'@importFrom parallel parLapply makeCluster clusterEvalQ clusterExport stopCluster
#'
#'@author Dustin Crockett
#'
#'@export
process_strong <- function(conn, core_count = 2, split_size=10, tolerance=0.053) {
  start_time <- Sys.time()
  
  
  #validate the conn list
  # verify that the conn is a list
  assertthat::assert_that(is.list(conn),
                          msg = "conn is not a list. Make sure conn is a list.")
  # verify that conn has the correct fields
  assertthat::assert_that("dbname" %in% names(conn),
                          msg = "dbname is missing from the list conn")
  assertthat::assert_that("host" %in% names(conn),
                          msg = "host is missing from the list conn")
  assertthat::assert_that("port" %in% names(conn),
                          msg = "port is missing from the list conn")
  assertthat::assert_that("user" %in% names(conn),
                          msg = "user is missing from the list conn")
  assertthat::assert_that("password" %in% names(conn),
                          msg = "password is missing from the list conn")
  
  con1 <- DBI::dbConnect(RMariaDB::MariaDB(),
                         dbname = conn$dbname,
                         host = conn$host,
                         port = conn$port,
                         user = conn$user,
                         password = conn$password)
  
  #calculate strong peptide
  results <- DBI::dbGetQuery(con1, "select id from candidate.genus")
  
  id_list <- split(results$id, rep(1:ceiling(length(results$id)/split_size), each=split_size)[1:length(results$id)])
  
  DBI::dbExecute(con1,'drop table if exists candidate.strong_peptide_temp')
  
  DBI::dbExecute(con1, 'create table candidate.strong_peptide_temp ( peptide_id int,
            genus_id int2,
            strong_percent float(23),
            primary key (peptide_id, genus_id) )')
  
  cl <- makeCluster(core_count)
  
  clusterEvalQ(cl, {
    ## set up each worker.  Could also use clusterExport()
    library(DBI)
    library(RMariaDB)
    NULL
  })
  #print(paste0("range ", range_test))
  full_result <- parLapply("cl"=cl, "X" = id_list, "fun" = analyze_strong_counts, chunk.size = 1, conn= conn)
  
  end_time = Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Populated strong_peptide_temp for all genus, started at ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds\n"))
  start_time = end_time
  
  stopCluster(cl)
  
  sqlcmd <- DBI::sqlInterpolate(con1,"insert into candidate.strong_peptides(peptide_id, strong_percent, strong_genus_id, second_percent, second_genus_id)
                           select id peptide_id, top_percent, top_genus_id, second_percent, second_genus_id  from
                           (select id, top_percent.strong_percent top_percent, top_percent.genus_id top_genus_id, second_percent.strong_percent second_percent, second_percent.genus_id second_genus_id
                           from candidate.peptides p
                           join candidate.strong_peptide_temp top_percent
                           on p.id = top_percent.peptide_id and top_percent.genus_id = (select genus_id from candidate.strong_peptide_temp where p.id = peptide_id order by strong_percent desc limit 1)
                           left join candidate.strong_peptide_temp second_percent
                           on p.id = second_percent.peptide_id and second_percent.genus_id = (select genus_id from candidate.strong_peptide_temp where p.id = peptide_id order by strong_percent desc limit 1,1)
                           ) results  where results.top_percent >= ifnull(results.second_percent,0) / ?tolerance or results.second_percent is null", tolerance = tolerance)
  DBI::dbExecute(con1, sqlcmd)
  
  
  end_time = Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Populated strong_peptides for all genus, started at ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds\n"))
  start_time = end_time
  
  #now we update the organisms table with the new strong_peptide_counts
  
  DBI::dbExecute(con1, "update candidate.organisms oo
	join (select o.id, count(distinct sp.peptide_id ) strong_pep_count
	from candidate.organisms o
		join candidate.species s on o.species_id = s.id
		join candidate.strong_peptides sp on sp.strong_genus_id = s.genus_id
		join candidate.organisms_to_peptides otp on otp.organism_id = o.id and otp.peptide_id = sp.peptide_id
	group by o.id) strong_counts on strong_counts.id = oo.id
	set oo.strong_peptide_count = strong_counts.strong_pep_count")
  
  #finally we drop the temp table
  #now that we no longer need the strong peptide temporary table we remove it
  DBI::dbExecute(con1, 'drop table candidate.strong_peptide_temp')
  
  DBI::dbDisconnect(con1)
  
  end_time = Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Populated strong_peptide_counts for all organisms, started at ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds\n"))
  
}
