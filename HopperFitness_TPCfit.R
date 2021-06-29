#Fit TPCs to M. dodgei hopping data

library(rTPC)

#hopping data
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/FitnessContrib_JEB/data/HopperTPCdata/")
jump.long= read.csv("HoppingData.csv")

#empirical TPC
#CONVERT FROM FT TO M
jump.long$dist =jump.long$dist*0.3048
#aggregate data
dat.pop= aggregate(jump.long, by=list(jump.long$Site, jump.long$Species, jump.long$temp), FUN=mean, na.action = na.omit)  #jump.long$Sex
names(dat.pop)[1:3]= c("Site","Species","temp")
dat.pop1= subset(dat.pop, dat.pop$Species=="dodgei" & dat.pop$elev==3048)

#add CTs to tpc data
dat.tpc= as.data.frame( cbind( c(dat.pop1[,"temp"],7.78,57.16), c(dat.pop1[,"dist"],0,0) ))
colnames(dat.tpc)=c("temp","rate") #change from dist
dat.tpc= dat.tpc[1:5,]

#fit loess
lo= loess(dist ~ temp, dat.tpc, span=1)

#---------------
#try rTPC package
#https://github.com/padpadpadpad/rTPC

nls_multstart(rate ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                         data = dat.tpc,
                         iter = 500,
                         start_lower = get_start_vals(dat.tpc$temp, dat.tpc$rate, model_name = 'lactin2_1995') - 2,
                         start_upper = get_start_vals(dat.tpc$temp, dat.tpc$rate, model_name = 'lactin2_1995') + 2,
                         supp_errors = 'Y')

dat.tpc$process="perf"
d_1=dat.tpc


d_models <- group_by(dat.tpc, process) %>%
  nest() %>%
  mutate(., thomas_2012 = map(data, ~nls_multstart(rate ~ thomas_2012(temp = temp, a, b, c, topt),
                                                data = .x,
                                                iter = 500,
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') - 1,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') + 2,
                                                supp_errors = 'Y',
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2012'))),
         gaussian = map(data, ~nls_multstart(rate ~ gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = 500,
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 2,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 2,
                                             supp_errors = 'Y')),
         quadratic = map(data, ~nls_multstart(rate ~ quadratic_2008(temp = temp, a, b, c),
                                              data = .x,
                                              iter = 500,
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 1,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 1,
                                              supp_errors = 'Y')),
         weibull = map(data, ~nls_multstart(rate ~ weibull_1995(temp = temp, a, topt, b, c),
                                            data = .x,
                                            iter = 1000,
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') -2,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') + 2,
                                            supp_errors = 'Y')),
         rezende = map(data, ~nls_multstart(rate ~ rezende_2019(temp = temp, a, q10, b, c),
                                            data = .x,
                                            iter = 500,
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 0.8,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 1.2,
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            supp_errors = 'Y')),
         beta = map(data, ~nls_multstart(rate ~ beta_2012(temp = temp, a, b, c, d, e),
                                         data = .x,
                                         iter = 500,
                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') -10,
                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') + 10,
                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         supp_errors = 'Y')),
         modgaussian = map(data, ~nls_multstart(rate ~ modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                                                data = .x,
                                                iter = 500,
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') - 1,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') +1,
                                                supp_errors = 'Y')),
         boatman = map(data, ~nls_multstart(rate ~ boatman_2017(temp = temp, rmax, tmin, tmax, a, b),
                                            data = .x,
                                            iter = 500,
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') -1,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 1,
                                            supp_errors = 'Y')),
         thomas_2017 = map(data, ~nls_multstart(rate ~ thomas_2017(temp = temp, a, b, c, d, e),
                                                data = .x,
                                                iter = 500,
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') *0.5,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') *1.5,
                                                supp_errors = 'Y')))


#===================
# stack models
d_stack <- gather(d_models, 'model', 'output', 6:ncol(d_models))

# preds
newdata <- tibble(temp = seq(min(d_1$temp), max(d_1$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(output, augment, newdata = newdata)) %>%
  unnest(preds)

# estimate parameters
params <- d_stack %>%
  mutate(., est = map(output, tidy)) %>%
  select(., -c(data, output)) %>%
  unnest(est)

# plot fit
ggplot(d_preds, aes(temp, rate)) +
  geom_point(aes(temp, rate), d_1) +
  geom_line(aes(temp, .fitted, col = model)) +
  facet_wrap(~model, labeller = labeller(model = label_facets_num)) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  xlab('Temperature (ºC)') +
  ylab('rate') +
  geom_hline(aes(yintercept = 0), linetype = 2)

#use rezende, function (temp, q10, a, b, c)
plot(10:40,rezende_2019(10:40, q10=2.27, a=0.109, b=9.02, c=0.00116), type="b")

nls_multstart(rate ~ rezende_2019(temp = temp, a, q10, b, c),
                                   data = dat.tpc,
                                   iter = 500,
                                   start_lower = get_start_vals(dat.tpc$temp, dat.tpc$rate, model_name = 'rezende_2019') * 0.8,
                                   start_upper = get_start_vals(dat.tpc$temp, dat.tpc$rate, model_name = 'rezende_2019') * 1.2,
                                   upper = get_upper_lims(dat.tpc$temp, dat.tpc$rate, model_name = 'rezende_2019'),
                                   lower = get_lower_lims(dat.tpc$temp, dat.tpc$rate, model_name = 'rezende_2019'),
                                   supp_errors = 'Y')
