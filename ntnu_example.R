table2005 <- to_table(c(4,2,20,0,0,7),2,3)
n_row2005 <- margins(table2005,1)
table2006 <- to_table(c(3,0,21,0,0,2),2,3)
n_row2006 <- margins(table2006,1)

table2007 <- to_table(c(4,0,23,0,2,5),2,3)
n_row2007 <- margins(table2007,1)

table2008 <- to_table(c(4,1,26,4,1,9),2,3)
n_row2008 <- margins(table2008,1)

table2009 <- to_table(c(0,0,40,5,2,11),2,3)
n_row2009 <- margins(table2009,1)

table2010 <- to_table(c(4,3,25,5,7,9),2,3)
n_row2010 <- margins(table2010,1)


fisher.test(table2005)

tables2005 <- gen_tables(n_row2005, 3)
reduced2005 <- group_reduce(tables2005)

ordering2005 <- supremum_ordering(n_row2005,
                                  3,
                                  show_progress = T, 
                                  N_order=15,
                                  pre_tables = tables2005,
                                  pre_group_reduced = reduced2005)
p2005 <- ordering2005[[2]][which(tables2005[[2]] == table_to_number(table2005))]


ordering2007 <- supremum_ordering(n_row2007,
                                  3,
                                  show_progress = T,
                                  N_order = 50,
                                  N_find = 100,
                                  pre_tables = tables2007,
                                  pre_group_reduced = reduced2007)