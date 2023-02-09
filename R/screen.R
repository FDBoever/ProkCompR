"histidine kinase"
"regulator"
prd <- "sensor"
prd <- "cytochrome"
prd <- "amylase"
prd <- "nitrate"
prd <- "lactate permease"
prd <- "permease"


prd <- "alginate"


annotation.tbl %>% filter(grepl(prd,product,ignore.case=T)) %>%
  select(genome,product) %>%
  group_by(genome) %>%
  tally() %>%
  arrange(desc(n)) %>%
  data.frame()


annotation.tbl %>% filter(grepl(prd,product,ignore.case=T)) %>%
  select(genome,product) %>%
  group_by(genome) %>%
  tally() %>%
  arrange(desc(n)) %>%
  data.frame()

annotation.tbl %>% filter(grepl(prd,product,ignore.case=T)) %>%
  select(genome,product) %>%
  group_by(genome,product) %>%
  tally() %>%
  arrange(desc(n)) %>%
  data.frame()

annotation.tbl %>% filter(grepl(prd,product,ignore.case=T)) %>%
  select(genome,product) %>%
  group_by(genome,product) %>%
  tally() %>%
  arrange(desc(genome)) %>%
  data.frame()
