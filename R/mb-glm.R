##### -- multiple general linear models (glm), broom and tidyverse

#Let's take a look at what we have
summary(df.chkm)

#select those that could be interesting...
#things like genomes, martkerrs and markerset could be ignored

df.glm = df.chkm %>%
  select(-c(BinId, Marker_lineage,genomes,markers,marker_sets, Translation_table, c0, c1, c2,c3,c4,c5plus))

library(broom)

lineReg <- glm(Genome_size ~.,data=df.glm)
summary(lineReg)


tidy(lineReg)
glance(lineReg)


aug.df.glm <- augment(lineReg, data=df.glm)


aug.df.glm %>% ggplot(aes(Genome_size,.fitted)) +
  geom_point(aes(color=group,size=.resid))

#--------------

#Or per group...
df.glm %>%
  #filter out the low-represented data? say Other????
  #filter(group!= 'other')
  group_by(group) %>%
  summarize(model = list(lm(Genome_size ~ Coding_density))) %>%
  mutate(tidied = map(model, tidy, conf.int = TRUE)) %>%
  unnest(tidied) %>%
  filter(term == 'Coding_density') %>%
  mutate(p.value.s = format.pval(p.value))



df.glm %>%
  group_by(group) %>%
  summarize(model = list(lm(Genome_size ~ Coding_density))) %>%
  mutate(tidied = map(model, tidy, conf.int = TRUE)) %>%
  unnest(tidied) %>%
  filter(term == 'Coding_density') %>%
  mutate(p.value.s = format.pval(p.value),
         group=fct_reorder(group, estimate)) %>%
  ggplot(aes(estimate, group)) +
    geom_point() +
    geom_errorbarh(aes(xmin=conf.low,xmax=conf.high),heigth = .1) +
    geom_vline(xintercept=0, lty=2) +
    labs(x='estimated slope',title='group specific relation coding density and genome size')





?map
#One could do this on a binary presence/absence matrix !
#binomial GLM, or get successes vs no successes type data...!

  #df.glm %>%
#  group_by(group) %>%

#  summarize(model = list(glm(Genome_size ~ Coding_density,family='binomial')))


