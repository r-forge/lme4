###################################################
### chunk number 1: preliminaries
###################################################
options(width=69,show.signif.stars=FALSE)
library(lattice)
lattice.options(default.theme = function() standard.theme())
#lattice.options(default.theme = function() standard.theme(color=FALSE))


###################################################
### chunk number 2: xyplotshow
###################################################
xyplot(optden ~ carb, Formaldehyde)


###################################################
### chunk number 3: xyplot
###################################################
print(xyplot(optden ~ carb, Formaldehyde, aspect = 1))


###################################################
### chunk number 4: xyplot1
###################################################
print(xyplot(optden ~ carb, Formaldehyde, aspect = 1,
             type = c("g","p","r"),
             xlab = "Amount of carbohydrate (ml)",
             ylab = "Optical density"))


###################################################
### chunk number 5: histshow eval=FALSE
###################################################
## histogram(~ count, InsectSprays)


###################################################
### chunk number 6: hist
###################################################
print(histogram(~ count, InsectSprays))


###################################################
### chunk number 7: densshow eval=FALSE
###################################################
## densityplot(~ count, InsectSprays)


###################################################
### chunk number 8: dens
###################################################
print(densityplot(~ count, InsectSprays))


###################################################
### chunk number 9: denssqrtshow eval=FALSE
###################################################
## densityplot(~ sqrt(count), InsectSprays,  xlab = "Square root of count")


###################################################
### chunk number 10: denssqrt
###################################################
print(densityplot(~ sqrt(count), InsectSprays,
                  xlab = "Square root of count"))


###################################################
### chunk number 11: denssqrtfancyshow eval=FALSE
###################################################
## densityplot(~ sqrt(count), InsectSprays, xlab = expression(sqrt("count")))


###################################################
### chunk number 12: denssqrtfancy
###################################################
print(densityplot(~ sqrt(count), InsectSprays,
                  xlab = expression(sqrt("count"))))


###################################################
### chunk number 13: comparedensshow eval=FALSE
###################################################
## densityplot(~ sqrt(count), InsectSprays, groups = spray,
##             auto.key = list(columns = 6))


###################################################
### chunk number 14: comparedens
###################################################
print(densityplot(~ sqrt(count), InsectSprays,
                  groups = spray,
                  auto.key = list(columns = 6),
                  xlab = expression(sqrt("count"))))


###################################################
### chunk number 15: comparedenspanelshow eval=FALSE
###################################################
## densityplot(~ sqrt(count)|spray, InsectSprays, layout = c(1,6))


###################################################
### chunk number 16: comparedenspanel
###################################################
print(densityplot(~ sqrt(count)|spray, InsectSprays,
                  layout = c(1,6),
                  xlab = expression(sqrt("count"))))


###################################################
### chunk number 17: comparedenspanelleftshow eval=FALSE
###################################################
## densityplot(~ sqrt(count)|spray, InsectSprays, layout = c(1,6), strip=FALSE,
##             strip.left = TRUE)


###################################################
### chunk number 18: comparedenspanelleft
###################################################
print(densityplot(~ sqrt(count)|spray, InsectSprays,
                  layout = c(1,6), strip=FALSE, strip.left = TRUE,
                  xlab = expression(sqrt("count"))))


###################################################
### chunk number 19: comparedenspanelreordershow eval=FALSE
###################################################
## densityplot(~ sqrt(count)|reorder(spray,count), InsectSprays)


###################################################
### chunk number 20: comparedensreorderpanel
###################################################
print(densityplot(~ sqrt(count)|reorder(spray, count), InsectSprays,
                  layout = c(1,6), strip=FALSE, strip.left = TRUE,
                  xlab = expression(sqrt("count"))))


###################################################
### chunk number 21: verticalbwshow eval=FALSE
###################################################
## bwplot(sqrt(count) ~ spray, InsectSprays)


###################################################
### chunk number 22: verticalbw
###################################################
print(bwplot(sqrt(count) ~ spray, InsectSprays,
             ylab = expression(sqrt("count"))))


###################################################
### chunk number 23: horbwshow eval=FALSE
###################################################
## bwplot(spray ~ sqrt(count), InsectSprays)


###################################################
### chunk number 24: horbw
###################################################
print(bwplot(spray ~ sqrt(count), InsectSprays,
             xlab = expression(sqrt("count"))))


###################################################
### chunk number 25: horbwreordshow eval=FALSE
###################################################
## bwplot(reorder(spray,count) ~ sqrt(count), InsectSprays)


###################################################
### chunk number 26: horreordbw
###################################################
print(bwplot(reorder(spray,count) ~ sqrt(count), InsectSprays,
             xlab = expression(sqrt("count"))))


###################################################
### chunk number 27: horbwcompshow eval=FALSE
###################################################
## bwplot(reorder(spray,count) ~ sqrt(count), InsectSprays, aspect = 0.2)


###################################################
### chunk number 28: horbwcomp
###################################################
print(bwplot(reorder(spray,count) ~ sqrt(count), InsectSprays, aspect = 0.2,
             xlab = expression(sqrt("count"))))


###################################################
### chunk number 29: dotsprayshow eval=FALSE
###################################################
## dotplot(reorder(spray,count) ~ sqrt(count), InsectSprays,
##         type = c("p","a"), pch = 21, jitter.y = TRUE)


###################################################
### chunk number 30: dotspray
###################################################
print(dotplot(reorder(spray,count) ~ sqrt(count), InsectSprays,
              aspect = 0.4, type = c("p","a"),
              pch = 21, jitter.y = TRUE,
              xlab = expression(sqrt("count"))))


