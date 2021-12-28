/* An efficient SAS macro for estimating direct adjusted survival functions
/* at the event times of K treatment groups for data with or without left truncation
/* based on a stratified Cox model, and for simulating pair-wise confidence bands
/* and simultaneous p-values.
/*
/* Author: Zhen-Huan Hu <novelkeny@gmail.com>
/*
/* Copyright (C) 2019-2021 Zhen-Huan Hu
/* -------------------------------------------- */
/* This program is free software: you can redistribute it and/or modify
/* it under the terms of the GNU General Public License as published by
/* the Free Software Foundation, either version 3 of the License, or
/* (at your option) any later version.
/*
/* This program is distributed in the hope that it will be useful,
/* but WITHOUT ANY WARRANTY; without even the implied warranty of
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
/* GNU General Public License for more details.
/*
/* You should have received a copy of the GNU General Public License
/* along with this program. If not, see <https://www.gnu.org/licenses/>.
/*
/* SAS and all other SAS Institute Inc. product or service names are registered trademarks
/* or trademarks of SAS Institute Inc. in the USA and other countries.
/* -------------------------------------------- */
/*
/* Function: %ADJSURVLT
/*
/* Parameter specifications:
/* -------------------------------------------- */
/*  INDATA:       Specifies the input data set, which requires the following variables:
/*                a continuous variable for the observed time, an event indicator, a categorical
/*                variable indicating treatment groups, and covariates for the multi-variate model.
/*                For left-truncated data, the input data set also should contain the time to
/*                left truncation. Observations with missing values need to be excluded beforehand.
/*  TIME:         Specifies the observed time from the start of follow-up to the event of interest
/*                or right censoring, whichever happens first.
/*  EVENT:        Specifies the event indicator variable.
/*                It takes 1 as event of interest and 0 as right censoring.
/*  STRATA:       Specifies the treatment group variable.
/*  COVLST:       Specifies the covariate list. Multiple covariates can be separated by spaces.
/*                A | sign can be used to separate continuous/ordinal variables from discrete
/*                categorical variables. For example, if variable A is categorical and variable B is
/*                continuous/ordinal, COVLST = A | B needs to be specified; if both are categorical,
/*                COVLST = A B needs to be specified; and if both are continuous/ordinal,
/*                COVLST = | A B needs to be specified.
/*  LTIME:        Specifies the left truncation time. The macro assumes no left truncation
/*                if it is left unspecified.
/*  SEED:         Specifies the seed for random number generation. By default, SEED = 0 which makes
/*                the seed obtained from the Intel RdRand instruction if feasible,
/*                or the current time value.
/*  NSIM:         Specifies the number of simulation processes for constructing confidence bands.
/*                By default, NSIM = 1000.
/*  MINTIME:      Specifies the lower boundary of the time interval between which the confidence bands
/*                are constructed. By default, the macro takes the minimum observed event time
/*                for each treatment group. For confidence bands of the survival difference estimates,
/*                the macro uses the greater of the minimum observed event times between the treatment
/*                group pairs. If the specified value is less than the minimum observed event time,
/*                the macro falls back to the default.
/*  MAXTIME:      Specifies the upper boundary of the time interval between which the confidence bands
/*                are constructed. By default, the macro takes the maximum observed event time
/*                for each treatment group. For confidence bands of the survival difference estimates,
/*                the macro uses the greater of the maximum observed event times between the treatment
/*                group pairs. If the specified value exceeds the maximum observed event time,
/*                the macro falls back to the default.
/*  ALPHA:        Specifies the confidence level. By default, ALPHA = 0.05.
/*  TIMELIST:     Specifies the time points at which to display the estimates.
/*                Setting it only affects the displayed results.
/*  OUTDATA:      Specifies the name of the main output data set containing direct adjusted survival
/*                estimates. It also defines the prefixes of the names of three additional output data
/*                sets that contain summary information of the sample, survival difference estimates
/*                and simultaneous p values.
/*                The default value is OUTDATA = ADJOUT, which creates the following data sets:
/*                [ADJOUT] which contains information associated with the adjusted survival estimates
/*                for individual treatment groups including event time, treatment groups,
/*                numbers at risk, survival estimates, standard error estimates, confidence limits
/*                and confidence bands;
/*                [ADJOUTINFO] which contains summary of the numbers of events and censored subjects
/*                as well as boundaries of the confidence bands;
/*                [ADJOUTDIFF] which contains information associated with the survival difference
/*                estimates between treatment group pairs including event time, treatment group pairs,
/*                survival difference estimates, standard error estimates, confidence limits
/*                and confidence bands;
/*                [ADJOUTTEST] which contains information regarding the hypothesis testing
/*                of the difference between treatment group pairs including number of simulation
/*                processes, critical values, and simultaneous p values.
/*                The confidence limits and confidence bands for the adjusted survival estimates are
/*                calculated based on log-log transformation, while those for the survival difference
/*                estimates are based on linear transformation.
/*  OUTSURVPLOT:  Enables the macro to plot direct adjusted survival curves.
/*                When OUTSURVPLOT = 1 is set, the macro outputs a plot file named ADJOUT-PLOT that
/*                displays a grid of adjusted survival curves, and when OUTSURVPLOT = 2 is set,
/*                a single plot with overlapped adjusted survival curves is produced instead.
/*  OUTDIFFPLOT:  Enables the macro to plot direct adjusted survival difference curves.
/*                When OUTDIFFPLOT = 1 is specified, an image file named ADJOUTDIFF-PLOT is created
/*                that displays a grid of survival difference curves.
/*  SHOWCI:       Toggles whether to display confidence limits in plots. The macro displays
/*                confidence limits by default. Setting SHOWCI = 0 hides the confidence limits.
/*  SHOWCB:       Toggles whether to display confidence bands in plots. The macro displays confidence
/*                bands by default. Setting SHOWCB = 0 hides the confidence bands.
/*  TICKVALUES:   Specifies values of the tick marks on the X axis.
/*  WIDTH:        Specifies width of the plots. By default, WIDTH = 800PX.
/*  HEIGHT:       Specifies height of the plots. By default, HEIGHT = 600PX.
/*  STYLE:        Specifies style of the plots. By default, STYLE = STATISTICAL.
/*  IMAGEFMT:     Specifies image format of the plots. By default, IMAGEFMT = PNG.
/*  NOPRINT:      Hides output results.
/* -------------------------------------------- */

%macro adjsurvlt(
  indata = , time = , event = , strata = , covlst = , ltime = ,
  seed = 0, nsim = 1000, mintime = , maxtime = , alpha = .05, timelist = ,
  outdata = adjout, outsurvplot = 0, outdiffplot = 0,
  showci = 1, showcb = 1, tickvalues = ,
  width = 800px, height = 600px, style = statistical, imagefmt = png,
  noprint = 0
  );

  %let dsid = %sysfunc(open(&indata, i));
  %let timenum = %sysfunc(varnum(&dsid, &time));
  %let timelbl = %qsysfunc(varlabel(&dsid, &timenum));
  %let eventnum = %sysfunc(varnum(&dsid, &event));
  %let eventlbl = %qsysfunc(varlabel(&dsid, &eventnum));
  %let stratanum = %sysfunc(varnum(&dsid, &strata));
  %let stratalbl = %qsysfunc(varlabel(&dsid, &stratanum));
  %let stratafmt = %qsysfunc(varfmt(&dsid, &stratanum));
  %let rc = %sysfunc(close(&dsid));
  %if &timelbl = %then %let timelbl = %upcase(&time);
  %if &eventlbl = %then %let eventlbl = %upcase(&event);
  %if &stratalbl = %then %let stratalbl = %upcase(&strata);
  %if &stratafmt = %then %let stratafmt = best8.;

  %if %sysfunc(prxmatch(/\|/, &covlst)) %then %do;
    %let catelst = %qsysfunc(prxchange(s/^(.*)\|(.*)$/$1/, 1, &covlst));
    %let contlst = %qsysfunc(prxchange(s/^(.*)\|(.*)$/$2/, 1, &covlst));
  %end;
  %else %do;
    %let catelst = &covlst;
    %let contlst = ;
  %end;

  proc iml;
    use &indata;
    %if &ltime = %then read all var {&time &event &strata} into xdata;
    %else read all var {&time &event &strata &ltime} into xdata;;
    %if &catelst ~= %then read all var {&catelst} into xcate;;
    %if &contlst ~= %then read all var {&contlst} into xcont;;
    close &indata;

    * Remove invalid values;
    miss_count = cmiss(xdata[, 3]);
    %if &catelst ~= %then miss_count = miss_count + (cmiss(xcate))[, +];;
    %if &contlst ~= %then miss_count = miss_count + (cmiss(xcont))[, +];;
    outidx = loc(miss_count = 0 & (xdata[, 2] = 1 | xdata[, 2] = 0) & xdata[, 1] > 0 %if &ltime ~= %then & xdata[, 4] >= 0;);

    xdata = xdata[outidx, ];
    %if &catelst ~= %then xcate = xcate[outidx, ];;
    %if &contlst ~= %then xcont = xcont[outidx, ];;

    * Create covariate matrix;
    nl = nrow(xdata);
    z = {};
    %if &catelst ~= %then %do;
      nx = ncol(xcate);
      do v = 1 to nx;
        var = xcate[, v];
        vdf = ncol(unique(var)) - 1; * DF;
        if vdf > 0 then do;
          z_design = design(var); * Encode levels in a design matrix;
          z = z || z_design[, 1: vdf];
        end;
      end;
    %end;
    %if &contlst ~= %then %do;
      z = z || xcont;
    %end;

    nz = ncol(z);
    if nz = 0 then do;
      z = j(nl, 1, 0);
      nz = 1;
    end;
    else if nz > 0 then do;
      call qr(q, r, piv, lindep, z); * QR decomposition;
      z = z[, piv[1: (nz - lindep)]]; * Remove linearly dependent columns;
      nz = ncol(z);
    end;
    znames = 'z' + strip(char(1: nz));
    store nz znames;

    out = xdata || z;
    %if &ltime = %then out_vnames = {'time' 'event' 'strata'} || znames;
    %else out_vnames = {'time' 'event' 'strata' 'ltime'} || znames;;
    create preproc from out[colname = out_vnames];
    append from out;
    close preproc;
  quit;

  proc phreg data = preproc covout outest = covout noprint;
    strata strata;
    %if &ltime = %then model time * event(0) = z:;
    %else model time * event(0) = z: / entry = ltime;;
    output out = coxout xbeta = zbeta;
  run;

  proc iml;
    load nz znames;
    use coxout;
    %if &ltime = %then read all var {time event strata zbeta};
    %else read all var {time event strata ltime zbeta};;
    read all var znames into z;
    close coxout;

    use covout where(_type_ = 'COV');
    read all var znames into covbeta; * Cov matrix for beta;
    close covout;

    nl = nrow(z); * N of z;
    byloc_event = loc(event);
    nevent = ncol(byloc_event); * N of events;
    bystrata = unique(strata);
    ns = ncol(bystrata); * N of strata;
    strata_ref = t(bystrata);
    strata_design = (strata = repeat(bystrata, nl, 1));

    * Get unique event times;
    pre_tlist = {0 &timelist &mintime &maxtime};
    tlist = t(unique(pre_tlist));
    ntlist = nrow(tlist);
    ts0 = j(ntlist # ns, 2, .);
    ts0[, 1] = repeat(tlist, ns, 1);
    ts0[, 2] = colvec(repeat(strata_ref, 1, ntlist));
    ts_block = ts0 // ((time || strata)[byloc_event, ]); * Event time only;

    call sortndx(tdx, ts_block, 2: 1);
    tdx_unique = uniqueby(ts_block, 2: 1, tdx);
    nt = nrow(tdx_unique); * N of unique event times;
    time_sorted = ts_block[tdx[tdx_unique], 1];
    strata_sorted = ts_block[tdx[tdx_unique], 2];
    strata_sorted_design = (strata_sorted = repeat(bystrata, nt, 1));
    free ts0 ts_block;

    * Get unique z;
    call sortndx(zdx, z, 1: nz);
    zdx_unique = uniqueby(z, 1: nz, zdx);
    nl_weighted = nrow(zdx_unique); * N of unique z;
    if nl_weighted > 1 then weight = (zdx_unique[2: nl_weighted] // (nl + 1)) - zdx_unique;
    else if nl_weighted = 1 then weight = nl;
    z_weighted = z[zdx[zdx_unique], ];
    zbeta_weighted = zbeta[zdx[zdx_unique]];

    * Create info matrix;
    infomat = j(ns, 3, .);
    do s = 1 to ns;
      sl = loc(strata_design[, s]);
      infomat[s, 2] = sum(event[sl] = 1); * N of event;
      infomat[s, 3] = sum(event[sl] = 0); * N of censored;
    end;
    infomat[, 1] = (infomat[, {2 3}])[, +];

    * Create risk matrix;
    s0beta = j(nevent, 1, .);
    s1beta = j(nevent, nz, .);
    do s = 1 to ns;
      sl = loc(strata_design[, s]);
      esl = loc(strata_design[byloc_event, s]);
      nsl = ncol(esl);
      %if &ltime = %then riskmat = ((time[byloc_event])[esl] <= repeat(t(time[sl]), nsl, 1));
      %else riskmat = ((time[byloc_event])[esl] <= repeat(t(time[sl]), nsl, 1) & (time[byloc_event])[esl] > repeat(t(ltime[sl]), nsl, 1));;
      s0beta[esl] = riskmat * exp(zbeta[sl]);
      s1beta[esl, ] = riskmat * (exp(zbeta[sl]) # z[sl, ]);
    end;
    ebeta = s1beta / s0beta; * t x E(beta);

    * Create integral matrix;
    atrisk = j(nt, 1, 0);
    histmat = j(nt, nevent, 0);
    do s = 1 to ns;
      sl = loc(strata_design[, s]);
      esl = loc(strata_design[byloc_event, s]);
      tl = loc(strata_sorted_design[, s]);
      ntl = ncol(tl);
      %if &ltime = %then atrisk[tl] = (time_sorted[tl] <= repeat(t(time[sl]), ntl, 1))[, +];
      %else atrisk[tl] = (time_sorted[tl] <= repeat(t(time[sl]), ntl, 1) & time_sorted[tl] > repeat(t(ltime[sl]), ntl, 1))[, +];;
      histmat[tl, esl] = (time_sorted[tl] >= repeat(t((time[byloc_event])[esl]), ntl, 1));
    end;

    * Calculate direct adjusted survivals;
    h0 = 1 / s0beta; * t x Baseline h(t);
    ch0 = histmat * h0; * t x Baseline cumulative h(t);
    s_weighted = exp(-(ch0 * t(exp(zbeta_weighted)))); * t x S(t|z);
    adjsurv = (s_weighted * weight) / nl; * t x Direct adjusted survivals;

    * Calculate SE;
    v1core = histmat * (h0 / s0beta);
    sszbeta = s_weighted * (weight # exp(zbeta_weighted));
    v1 = (sszbeta ## 2) # v1core / (nl ## 2);

    if nl_weighted > 1 then do;
      v2core = histmat * (ebeta # h0);
      zzbeta_weighted = z_weighted # exp(zbeta_weighted);
      ssh = j(nt, nz, .); * t x Sum of S(t|z) multiplied by H matrix;
      do c = 1 to nz;
        ssh[, c] = (s_weighted # (
        ch0 * t(zzbeta_weighted[, c]) - v2core[, c] * t(exp(zbeta_weighted))
        )) * weight;
      end;
      sshcovbeta = ssh * covbeta;
      v2 = (sshcovbeta # ssh)[, +] / (nl ## 2);
      adjvar = v1 + v2; * t x Variances;
    end;
    else if nl_weighted = 1 then adjvar = v1;
    adjse = sqrt(adjvar); * t x SE;
    free s_weighted zzbeta_weighted zbeta_weighted z_weighted;

    * Calculate CL - log-log transformation;
    z_alpha = probit(1 - &alpha / 2.0);
    cl = j(nt, 2, .);
    byloc_ci = loc(adjsurv > 0 & adjsurv < 1);
    if col(byloc_ci) > 0 then do;
      tor = adjse[byloc_ci] / (adjsurv[byloc_ci] # log(adjsurv[byloc_ci]));
      cl[byloc_ci, 1] = adjsurv[byloc_ci] ## exp(-z_alpha # tor);
      cl[byloc_ci, 2] = adjsurv[byloc_ci] ## exp(z_alpha # tor);
    end;

    * Construct time mapping;
    sc = allcomb(ns, 2);
    call sort(sc, 1: 2);
    nsc = nrow(sc);
    strata_pair_ref = shape(bystrata[shape(sc, 1)], nsc, 2);

    time_pooled = {};
    sc_pooled = {};
    byloc = {}; * Mapping vectors;
    byloc_offset = loc(time_sorted = 0) - 1;
    do s = 1 to nsc;
      time_sorted_s1 = time_sorted[loc(strata_sorted_design[, sc[s, 1]])];
      time_sorted_s2 = time_sorted[loc(strata_sorted_design[, sc[s, 2]])];
      time_pooled_s12 = t(unique(time_sorted[loc((strata_sorted_design[, sc[s, ]])[, +] &
      time_sorted <= max(time_sorted_s1, time_sorted_s2))]));
      time_pooled = time_pooled // time_pooled_s12;

      nt_pooled_s12 = nrow(time_pooled_s12);
      sc_pooled = sc_pooled // j(nt_pooled_s12, 1, s);

      byloc_s12 = j(nt_pooled_s12, 2, .);
      byloc_s12[, 1] = (repeat(t(time_sorted_s1), nt_pooled_s12, 1) <= time_pooled_s12)[, +] + byloc_offset[sc[s, 1]];
      byloc_s12[, 2] = (repeat(t(time_sorted_s2), nt_pooled_s12, 1) <= time_pooled_s12)[, +] + byloc_offset[sc[s, 2]];
      byloc = byloc // byloc_s12;
    end;
    nt_pooled = nrow(time_pooled);
    strata_pair = j(nt_pooled, 2, .);
    strata_pair[, 1] = strata_sorted[byloc[, 1]];
    strata_pair[, 2] = strata_sorted[byloc[, 2]];

    * Calculate S1(t) - S2(t) and its SE;
    wdiff = adjsurv[byloc[, 1]] - adjsurv[byloc[, 2]];
    if nl_weighted > 1 then do;
      thetacovbeta = sshcovbeta[byloc[, 1], ] - sshcovbeta[byloc[, 2], ];
      theta = ssh[byloc[, 1], ] - ssh[byloc[, 2], ];
      vardiff = v1[byloc[, 1]] + v1[byloc[, 2]] + (thetacovbeta # theta)[, +] / (nl ## 2);
    end;
    else if nl_weighted = 1 then vardiff = v1[byloc[, 1]] + v1[byloc[, 2]];
    sediff = sqrt(vardiff);
    free thetacovbeta theta;

    * Calculate CL for S1(t) - S2(t);
    cldiff = wdiff + (z_alpha # sediff) * {-1 1};

    * Initialize simulations for CB;
    nsim = &nsim;
    gmat = j(nl, nsim, .);
    call randseed(&seed);
    call randgen(gmat, 'NORMAL');
    gevent = gmat[byloc_event, ];

    w1_sim = -(sszbeta # (histmat * (gevent / s0beta))) / nl;
    if nl_weighted > 1 then do;
      w2_sim = -(sshcovbeta * ((t(z[byloc_event, ] - ebeta)) * gevent)) / nl;
      w_sim_add = w1_sim + w2_sim;
      w_sim_sub = w1_sim - w2_sim;
    end;
    else if nl_weighted = 1 then do;
      w_sim_add = w1_sim;
      w_sim_sub = w1_sim;
    end;
    free histmat gevent gmat w1_sim w2_sim;

    * Define t1, t2 for CB;
    t1_prebyloc = j(nt, 1, .);
    t2_prebyloc = j(nt, 1, .);
    do s = 1 to ns;
      min_et = min(time[loc(strata_design[, s] & event = 1)]); * Min event time;
      max_et = max(time[loc(strata_design[, s] & event = 1)]); * Max event time;
      tl = loc(strata_sorted_design[, s]);
      t1_prebyloc[tl] = (time_sorted[tl] >= min_et);
      t2_prebyloc[tl] = (time_sorted[tl] <= max_et);
    end;
    t1diff_prebyloc = t1_prebyloc[byloc[, 1]] & t1_prebyloc[byloc[, 2]];
    t2diff_prebyloc = t2_prebyloc[byloc[, 1]] | t2_prebyloc[byloc[, 2]];

    %if &mintime ~= %then %do;
      t1_prebyloc = t1_prebyloc & (time_sorted >= &mintime);
      t1diff_prebyloc = t1diff_prebyloc & (time_pooled >= &mintime);
    %end;
    %if &maxtime ~= %then %do;
      t2_prebyloc = t2_prebyloc & (time_sorted <= &maxtime);
      t2diff_prebyloc = t2diff_prebyloc & (time_pooled <= &maxtime);
    %end;

    t12_prebyloc = t1_prebyloc & t2_prebyloc;
    t12diff_prebyloc = t1diff_prebyloc & t2diff_prebyloc;

    * Calculate CB for S(t) - log-log transformation;
    cb = j(nt, 2, .);
    t12mat = j(ns, 2, .);
    do s = 1 to ns;
      t12_byloc = loc(strata_sorted_design[, s] & t12_prebyloc);
      if col(t12_byloc) > 0 then do;
        t12mat[s, 1] = min(time_sorted[t12_byloc]);
        t12mat[s, 2] = max(time_sorted[t12_byloc]);

        qb = t((abs(w_sim_add[t12_byloc, ] / adjse[t12_byloc]))[<>, ]);
        call qntl(c_alpha, qb, 1 - &alpha, 2);
        torcb = adjse[t12_byloc] / (adjsurv[t12_byloc] # log(adjsurv[t12_byloc]));
        cb[t12_byloc, 1] = adjsurv[t12_byloc] ## exp(-c_alpha # torcb);
        cb[t12_byloc, 2] = adjsurv[t12_byloc] ## exp(c_alpha # torcb);
      end;
    end;

    * Calculate CB for S1(t) - S2(t);
    cp = j(nsc, 3, .);
    cp[, 1] = nsim;

    cbdiff = j(nt_pooled, 2, .);
    wdiff_sim = w_sim_sub[byloc[, 1], ] + w_sim_add[byloc[, 2], ];
    do s = 1 to nsc;
      t12diff_byloc = loc(sc_pooled = s & t12diff_prebyloc);
      if col(t12diff_byloc) > 0 then do;
        qb = t((abs(wdiff_sim[t12diff_byloc, ] / sediff[t12diff_byloc]))[<>, ]);
        call qntl(c_alpha, qb, 1 - &alpha, 2); * Critical value;
        cbdiff[t12diff_byloc, ] = wdiff[t12diff_byloc] + (c_alpha # sediff[t12diff_byloc]) * {-1 1};

        cp[s, 2] = c_alpha;
        qmax = max(abs(wdiff[t12diff_byloc] / sediff[t12diff_byloc]));
        cp[s, 3] = mean(qb > qmax); * P-value;
      end;
    end;
    free w_sim_add w_sim_sub wdiff_sim;

    outinfo = strata_ref || infomat || t12mat;
    create &outdata.info from outinfo[colname = {"&strata" 'total' 'event' 'censored' 't1' 't2'}];
    append from outinfo;
    close &outdata.info;

    outsurv = time_sorted || strata_sorted || adjsurv || adjse || cl || cb || atrisk;
    create &outdata from outsurv[colname = {"&time" "&strata" 'adjsurv' 'stderr' 'lcl' 'ucl' 'lcb' 'ucb' 'atrisk'}];
    append from outsurv;
    close &outdata;

    outdiff = time_pooled || strata_pair || wdiff || sediff || cldiff || cbdiff;
    create &outdata.diff from outdiff[colname = {"&time" "&strata.1" "&strata.2" 'adjdiff' 'stderr' 'lcl' 'ucl' 'lcb' 'ucb'}];
    append from outdiff;
    close &outdata.diff;

    outtest = strata_pair_ref || cp;
    create &outdata.test from outtest[colname = {"&strata.1" "&strata.2" 'nsim' 'c_alpha' 'p_value'}];
    append from outtest;
    close &outdata.test;
  quit;

  %if &noprint = 0 %then %do;
    %if &timelist = %then %let whereopt = ;
    %else %let whereopt = (where = (&time in (&timelist)));
    proc report data = &outdata.info nowd split = '|' headline;
      title "Summary of Event and Censoring";
      column &strata total event censored ("Band Boundaries" t1 t2);
      define &strata / display "&stratalbl" format = &stratafmt;
      define total / display "Total" format = 6.0;
      define event / display "Event" format = 6.0;
      define censored / display "Censor" format = 6.0;
      define t1 / display "Lower" format = best6.2;
      define t2 / display "Upper" format = best6.2;
    run;
    proc report data = &outdata &whereopt nowd split = '|' headline;
      title "Direct Adjusted Survival Estimates";
      column &time &strata adjsurv stderr ("Conf Limits" lcl ucl) ("Conf Bands" lcb ucb) atrisk;
      define &time / display "&timelbl" format = best6.2;
      define &strata / display "&stratalbl" format = &stratafmt;
      define adjsurv / display "Surv" format = 6.4;
      define stderr / display "Std Err" format = 6.4;
      define lcl / display "Lower" format = 6.4; define ucl / display "Upper" format = 6.4;
      define lcb / display "Lower" format = 6.4; define ucb / display "Upper" format = 6.4;
      define atrisk / display "N at Risk" format = 6.0;
    run;
    proc report data = &outdata.diff &whereopt nowd split = '|' headline;
      title "Direct Adjusted Survival Difference Estimates";
      column &time &strata.1 &strata.2 adjdiff stderr ("Conf Limits" lcl ucl) ("Conf Bands" lcb ucb);
      define &time / display "&timelbl" format = best6.2;
      define &strata.1 / display "&stratalbl 1" format = &stratafmt;
      define &strata.2 / display "&stratalbl 2" format = &stratafmt;
      define adjdiff / display "Surv Diff" format = 6.4;
      define stderr / display "Std Err" format = 6.4;
      define lcl / display "Lower" format = 6.4; define ucl / display "Upper" format = 6.4;
      define lcb / display "Lower" format = 6.4; define ucb / display "Upper" format = 6.4;
    run;
    proc report data = &outdata.test nowd split = '|' headline;
      title "Hypothesis Tests";
      column &strata.1 &strata.2 nsim c_alpha p_value;
      define &strata.1 / display "&stratalbl 1" format = &stratafmt;
      define &strata.2 / display "&stratalbl 2" format = &stratafmt;
      define nsim / display "N of Sim" format = 6.0;
      define c_alpha / display "Crit Value" format = 6.4;
      define p_value / display "P Value" format = pvalue6.4;
    run;
  %end;

  %if &outsurvplot > 0 or &outdiffplot > 0 %then %do;
    %let conf = %sysevalf(100 * (1 - &alpha));
    %if &mintime = %then %let minopt = ; %else %let minopt = min = &mintime;
    %if &maxtime = %then %let maxopt = ; %else %let maxopt = max = &maxtime;
    %if &tickvalues = %then %let valuesopt = ;
    %else %let valuesopt = values = (&tickvalues);
    %if &outsurvplot > 0 %then %do;
      ods listing style = &style;
      ods graphics on / imagename = "&outdata-plot" imagefmt = &imagefmt width = &width height = &height;
      title "Direct Adjusted Survival Curves for &eventlbl";
      %if &outsurvplot = 1 %then %do;
        proc sgpanel data = &outdata;
          panelby &strata / skipemptycells;
          %if &showcb %then band x = &time lower = lcb upper = ucb / name = 'cb' type = step fill transparency = .5 legendlabel = "&conf% Simultaneous CB";;
          %if &showci %then band x = &time lower = lcl upper = ucl / name = 'cl' type = step nofill legendlabel = "&conf% Pointwise CI";;
          step x = &time y = adjsurv / name = 'adjsurv' justify = left legendlabel = "Survival Estimates";
          colaxis label = "&timelbl" &minopt &maxopt &valuesopt;
          rowaxis label = "Direct Adjusted Survival Probability" min = 0 max = 1;
          keylegend 'adjsurv' 'cl' 'cb' / position = bottom;
          label &strata = "&stratalbl";
          format &strata &stratafmt;
        run;
      %end;
      %else %if &outsurvplot = 2 %then %do;
        proc sgplot data = &outdata;
          %if &showcb %then band x = &time lower = lcb upper = ucb / name = 'cb' group = &strata type = step fill transparency = .5 legendlabel = "&conf% Simultaneous CB";;
          %if &showci %then band x = &time lower = lcl upper = ucl / name = 'cl' group = &strata type = step nofill legendlabel = "&conf% Pointwise CI";;
          step x = &time y = adjsurv / name = 'adjsurv' group = &strata justify = left;
          xaxis label = "&timelbl" &minopt &maxopt &valuesopt;
          yaxis label = "Direct Adjusted Survival Probability" min = 0 max = 1;
          keylegend 'adjsurv' 'cl' 'cb' / position = bottom;
          label &strata = "&stratalbl";
          format &strata &stratafmt;
        run;
      %end;
      ods graphics off;
    %end;
    %if &outdiffplot = 1 %then %do;
      proc sort data = &outdata.test; by &strata.1 &strata.2;
      proc sort data = &outdata.diff; by &strata.1 &strata.2;
      data outdiffplot;
        merge &outdata.test &outdata.diff;
        by &strata.1 &strata.2;
        length panellabel $ 250;
        panellabel = catx(': ', "&stratalbl", catx(' - ', put(&strata.1, &stratafmt), put(&strata.2, &stratafmt)));
      run;
      ods listing style = &style;
      ods graphics on / imagename = "&outdata.diff-plot" imagefmt = &imagefmt width = &width height = &height;
      title "Differences between Direct Adjusted Survival Curves for &eventlbl";
      proc sgpanel data = outdiffplot;
        panelby panellabel / novarname skipemptycells;
        %if &showcb %then band x = &time lower = lcb upper = ucb / name = 'cb' type = step fill transparency = .5 legendlabel = "&conf% Simultaneous CB";;
        %if &showci %then band x = &time lower = lcl upper = ucl / name = 'cl' type = step nofill legendlabel = "&conf% Pointwise CI";;
        step x = &time y = adjdiff / name = 'adjdiff' justify = left legendlabel = "Survival Difference Estimates";
        colaxis label = "&timelbl" &minopt &maxopt &valuesopt;
        rowaxis label = "Difference of Survival Probabilities" offsetmin = 0.1 offsetmax = 0.1;
        refline 0 / axis = y;
        inset p_value / position = topleft;
        keylegend 'adjdiff' 'cl' 'cb' / position = bottom;
        label p_value = "Simultaneous P Value";
        format p_value pvalue6.4;
      run;
      ods graphics off;
    %end;
  %end;
%mend;
