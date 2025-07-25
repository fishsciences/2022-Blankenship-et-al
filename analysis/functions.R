# Data cleaning and convenience functions used in analysis
# M. Johnston, M. Espe
# Tue Sep 28 14:45:18 2021 ------------------------------

# quickly 'vet' a dataframe by previewing rows at the head, middle, and tail:
#-------------------------------------------------------#
vet <- function(d, n = 4L) {
 if(class(d) != 'data.frame') stop('vet() can only vet dataframes')
 left <- as.integer(nrow(d) / 2 - n / 2)
 torso = d[seq_len(n) + left - 1L,]
rbind(head(d, n), torso, tail(d, n))
}


# typing shortcuts
#--------------------------------------------#
len <- function(x){length(unique(x))}
csn <- function(x){colSums(is.na(x))}
rsn <- function(x){rowSums(is.na(x))}


# Hydro data functions
#-------------------------------------------------------#
mk_treatments = function(df)
{
    cc = colnames(df)
    cc = cc[grep("NEG|POS", cc)]
    distances = as.integer(gsub("^[A-Z]{3}", "", cc))

    side = factor(gsub("[0-9]+$", "", cc))
    distances[side == "NEG"] = distances[side == "NEG"] * -1
    levels(side) = c("upstream","downstream")
    return(data.frame(Side = side,
                      Distance_m = distances,
                      column = match(cc, colnames(df))))
}

hydro_vals = function(tidal, hydro, time_col = "DateTimePT",
                      merge_1m = TRUE)
{
    trt = mk_treatments(hydro)
    if(merge_1m){
        trt = rbind(trt, c(Side = NA,
                           Distance_m = 0,
                           column = grep("Car", colnames(hydro))))
        tidal$Distance_m[abs(tidal$Distance_m) == 1] = 0
    }

    rows = match(tidal[[time_col]], hydro$Time)
    cols = seq(nrow(trt))[match(tidal$Distance_m, trt$Distance_m)]
    as.matrix(hydro[,trt$column])[cbind(rows, cols)]
}

