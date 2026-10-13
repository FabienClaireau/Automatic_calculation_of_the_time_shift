#####################################################
###        Acoustic Trajectory - AFPR Approach    ###
###          Script 2 - new decalage auto         ###
###                     ######                    ###
### Script F. Claireau & Y. Bas - New formatting  ###
###           & corrections T. Ferraille          ###
#####################################################

####################
### Introduction ###
####################

# The acoustic trajectory approach using AFPR, developed by Fabien CLAIREAU and Yves BAS of CESCO-MNHN: Claireau F. et al (2018) Bat overpasses: an insufficient solution to encourage bats across roads https://doi.org/10.1111/1365-2664.13288 MNHN-CESCO/Greifswald University/Naturalia Environnement

# Three R scripts must be run successively:
# -   Script 1 - New Pairing: calculates the TDOA between two microphones of a recorder
# -   Script 2 - New Auto-Offset: calculates the time offset between two recorders
# -   Script 3 - New Find Bat Crossings: matches two recorders and defines crossings

## Inputs and Outputs for Script 2:

# **Inputs for Script 2:**
# -   ListPairs.csv: Fill out the file and save it as a CSV file (delimiter: semicolon); enter Inter_distance in meters.
# -   TrajTot.csv: Output from Script 1.

# **Outputs for Script 2:**
# -   EXPORTED FILE: ListPairs_Dec.csv
# -   Open the file and check if the 'decalageFinal' column contains any NA values, which would indicate that one of the recorders in the pair recorded in mono rather than stereo.
# -   Delete the row if this is the case.
# -   In any case, save the file under the name 'ListPairs_Dec_Man' in CSV format (delimiter: semicolon).



####################################################
### Parameter Section to Modify: paths and files ###
####################################################

rm(list = ls()) # Clear existing objects in the workspace

# Set the path to your tables (excel files, etc (WARNING use "/" instead of "\")
chemin_travail <- "C:/Users/FabienCLAIREAU/Desktop/Formation_MNHN/tables"



########################################
### Execution Section: DO NOT MODIFY ###
###       CTRL A > CTRL ENTER        ###
########################################


# Installing and loading packages
install_load_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}
install_load_packages(c("data.table", "Hmisc","dplyr"))

# Loading CSV files
setwd(chemin_travail)
TrajTot <- read.csv2("TrajTot.csv")
ListPaires <- read.csv2("ListPairs.csv")

# Initialisation of the output variable
decalageFinal <- vector()

# Maximum offset value
maxdecalage <- 5  

# === Calculation of decages per pair of participants === #

for (h in 1:nrow(ListPaires)) {
  
  # Filter the data for the first participant in the pair
  TrajTot1 <- subset(TrajTot, participation %in% ListPaires$participation_1[h])
  if (nrow(TrajTot1) == 0) {
    decalageTemp <- NA
  } else {
    TrajTot1 <- cbind(TrajTot1, Cote = 1, Point = ListPaires$point_1[h])
    
    # Filter the data for the second participant in the pair
    TrajTot2 <- subset(TrajTot, participation %in% ListPaires$participation_2[h])
    if (nrow(TrajTot2) == 0) {
      decalageTemp <- NA
    } else {
      TrajTot2 <- cbind(TrajTot2, Cote = 2, Point = ListPaires$point_2[h])
      
      # Merge the two participants for comparison
      Trav <- rbind(TrajTot1, TrajTot2)
      
      # Converting columns into numeric format
      Trav$DecDeb <- as.numeric(as.character(Trav$DecDeb))
      Trav$DecFin <- as.numeric(as.character(Trav$DecFin))
      Trav$Prob <- as.numeric(as.character(Trav$probabilite))
      
      # Extracting time information from file names
      Date <- vector(length = nrow(Trav))
      Heure <- vector(length = nrow(Trav))
      Minute <- vector(length = nrow(Trav))
      Seconde <- vector(length = nrow(Trav))
      MiliSec <- vector(length = nrow(Trav))
      Fich <- as.character(Trav$File0)
      
      for (i in 1:length(Fich)) {
        Date[i] <- substr(Fich[i], nchar(Fich[i]) - 18, nchar(Fich[i]) - 11)
        Heure[i] <- substr(Fich[i], nchar(Fich[i]) - 9, nchar(Fich[i]) - 8)
        Minute[i] <- substr(Fich[i], nchar(Fich[i]) - 7, nchar(Fich[i]) - 6)
        Seconde[i] <- substr(Fich[i], nchar(Fich[i]) - 5, nchar(Fich[i]) - 4)
        MiliSec[i] <- substr(Fich[i], nchar(Fich[i]) - 2, nchar(Fich[i]))
      }
      
      # Converting temporal elements into digital
      Datenum <- as.numeric(Date)
      Heure <- as.numeric(Heure)
      Minute <- as.numeric(Minute)
      Seconde <- as.numeric(Seconde)
      MiliSec <- as.numeric(MiliSec)
      Temps <- Heure * 3600 + Minute * 60 + Seconde + MiliSec / 1000
      
      # Adding time columns to the table
      Trav <- cbind(Trav, Heure, Minute, Seconde, MiliSec, Temps, Date, Datenum)
      
      # Selection of inputs/outputs for each participant
      Ent1 <- subset(Trav, DecDeb > 0 & DecFin < 0 & Cote == 1)
      Sort1 <- subset(Trav, DecDeb < 0 & DecFin > 0 & Cote == 1)
      Ent2 <- subset(Trav, DecDeb > 0 & DecFin < 0 & Cote == 2)
      Sort2 <- subset(Trav, DecDeb < 0 & DecFin > 0 & Cote == 2)
      
      # Calculating correspondences between participants
      if ((((nrow(Ent1) == 0) | nrow(Sort2) == 0)) & ((nrow(Ent2) == 0) | (nrow(Sort1) == 0))) {
        decalageTemp <- 0
      } else {
        Interval <- 5
        NbMatch <- vector()
        
        # Loop on offset intervals
        for (j in (-(round(maxdecalage / Interval)):round(maxdecalage / Interval))) {
          decalage <- j * Interval
          EntTemp <- Ent1$Temps + decalage
          SortTemp <- Sort1$Temps + decalage
          NbMatchTemp <- 0
          
          # Correspondence for inputs and outputs
          if (((nrow(Ent1) > 0) & (nrow(Sort2) > 0))) {
            L1 <- levels(as.factor(Ent1$espece))
            L2 <- levels(as.factor(Sort2$espece))
            L12 <- unique(c(L1, L2))
            SpNum1 <- match(Ent1$espece, L12)
            SpNum2 <- match(Sort2$espece, L12)
            TC12 <- find.matches(as.matrix(cbind(Ent1$Datenum, SpNum1, EntTemp)),
                                 as.matrix(cbind(Sort2$Datenum, SpNum2, Sort2$Temps)),
                                 tol = c(0, 0, (ListPaires$inter_distance[h] / 5 + (ListPaires$inter_distance[h] / 5 * 0.5) + 5)),
                                 maxmatch = 1)
            
            NbMatchTemp <- NbMatchTemp + length(subset(TC12$matches, TC12$matches > 0))
          }
          if (((nrow(Ent2) > 0) & (nrow(Sort1) > 0))) {
            L1 <- levels(as.factor(Ent2$espece))
            L2 <- levels(as.factor(Sort1$espece))
            L12 <- unique(c(L1, L2))
            SpNum1 <- match(Ent2$espece, L12)
            SpNum2 <- match(Sort1$espece, L12)
            TC21 <- find.matches(cbind(Ent2$Datenum, SpNum1, Ent2$Temps),
                                 cbind(Sort1$Datenum, SpNum2, SortTemp),
                                 tol = c(0, 0, (ListPaires$inter_distance[h] / 5 + (ListPaires$inter_distance[h] / 5 * 0.5) + 5)),
                                 maxmatch = 1)
            NbMatchTemp <- NbMatchTemp + length(subset(TC21$matches, TC21$matches > 0))
          }
          
          # Store the number of matches found
          NbMatch <- c(NbMatch, NbMatchTemp)
          print(paste(h, "/", nrow(ListPaires), " ", j * Interval))
        }
        
        # Selection of offsets with the most matches
        DecPotentiels <- (-(round(maxdecalage / Interval)):round(maxdecalage / Interval)) * Interval
        barplot(NbMatch, names.arg = DecPotentiels, cex.names = 0.5, main = h)
        DecSelectionnes <- subset(DecPotentiels, NbMatch == max(NbMatch))
        decalageTemp <- quantile(DecSelectionnes, 0.5, type = 1)
      }
    }
  }
  
  # Storage of the final offset for each pair
  decalageFinal <- c(decalageFinal, decalageTemp)
}

# === Saving results === #
ListDecalage <- cbind(ListPaires, decalageFinal)
write.csv2(ListDecalage, "ListPairs_Dec.csv")