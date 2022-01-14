
stmf <- read.csv("https://www.mortality.org/Public/STMF/Outputs/stmf.csv",header=TRUE,sep=",",skip=2)
stmf$CountryCode <- as.factor(stmf$CountryCode)
levels(stmf$CountryCode)[1] <- "AUS"

#Function to read HMD life table
hmd.lt <- function (country, username, password, label = country) 
{
    path <- paste("https://www.mortality.org/hmd/", country, 
        "/STATS/", "bltper_1x1.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    life.table <- try(utils::read.table(con, skip = 2, header = TRUE, 
        na.strings = "."), TRUE)
    close(con)
    if (class(life.table) == "try-error") 
        stop("Connection error at www.mortality.org. Please check username, password and country label.")
    return(life.table)
}

nax <- function(lt,age){  # determine nax from HMD 1x1 LT for age interval [x1,xn]
    # Credits to Carl Boe
  j1=age[-length(age)]+1; jn=age[-1]+1;
  nLx =( lt$T[j1] - lt$T[jn]);
  ndx=lt$l[j1]-lt$l[jn] ;
  nax = (nLx - (jn-j1)* lt$l[jn])/ ndx;
  return(nax)
  }

#Simple function to get life table from rates,
get.lt <- function(mx,age,sex="F",ax){ 
    ageint <- age[-1]-age[-length(age)]
    ax <- as.numeric(c(ax,ifelse(mx[length(age)]>0,1/mx[length(age)],0)))
    qx <- rep(NA,length(age))
    qx[1:length(age)-1] <- as.numeric(ageint*mx[1:length(age)-1]/(1+(ageint-ax[1:length(age)-1])*mx[1:length(age)-1]))
    qx[length(age)] <- 1
    lx <- rep(NA,length(age))
    lx <- round(c(100000,100000*cumprod(1-qx[-length(age)])))
    lx[lx<0|is.na(lx)] <- 0
    dx <- c(-diff(lx),lx[length(lx)])
    Lx <- (c(ageint,110-age[length(age)]))*lx-(c(ageint,110-age[length(age)])-ax)*dx
    Lx[length(age)] <- sum(c(0,lx[length(age)]*ax[length(age)]),na.rm=TRUE)
    Tx <- rev(cumsum(rev(Lx)))
    ex <- ifelse(lx > 0, Tx/lx , NA)
    
    return(data.frame(x=age,nax=ax,nmx=mx,nqx=qx,lx=lx,
                     ndx=dx,nLx=Lx,Tx=Tx,ex=ex))
}

age <- c(0,15,65,75,85)

mymatrix <- data.frame("Country"=NA,"Year"=NA,"week"=NA,"e0.week"=NA,"e0.hmd"=NA)

for (i in 2:length(names(table(stmf$CountryCode)))){
    
    cnt <- names(table(stmf$CountryCode))[i]

    #get input data, to include deaths with missing dates 
    input <- read.csv(paste("STMFinput/",cnt,"stmf.csv",sep=""),header=TRUE,sep=",",na.strings=".")
        
        #GET life table from HMD
    cnt.lt <- hmd.lt(country=cnt,user="stefano.mazzuco@unipd.it",password="08122018")
    range.years <- range(as.numeric(names((table(stmf$Year[stmf$CountryCode == cnt])))))

    for(y in range.years[1]:range.years[2]){
        D.y <- sum(stmf[which(stmf$CountryCode==cnt&stmf$Sex=="b"&stmf$Year==y),10],na.rm=TRUE)
        D.nodate <- input$Deaths[which(input$Sex=="b"&input$Year==y&input$Week=="UNK"&input$Age=="TOT")]
        D.nodate <- ifelse(length(D.nodate)==0,0,D.nodate)
        Av.Weeks <- max(as.numeric(names(table(stmf$Week[stmf$Country==cnt&stmf$Year==y]))))

        min.Weeks <- min(as.numeric(names(table(stmf$Week[stmf$Country==cnt&stmf$Year==y]))))
        for(w in min.Weeks:Av.Weeks){
            ##RATES. Note I adjust rates by including deaths with missing dates. I'm assuming that deaths with mssing dates are uniformly distribued across age classes,
            ##There are some countries (e.g. Sweden) where the number of deaths with missing dates is relatively high
            mx.y <- ((D.nodate+D.y)/D.y)*stmf[which(stmf$CountryCode==cnt&stmf$Sex=="b"&stmf$Year==y&stmf$Week==w),11:15]
            m.hmd <- c(sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 0:14)])/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 0:14)]),
                sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 15:64)])/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 15:64)]),
                sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 65:74)])/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 65:74)]),
                sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 75:84)])/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 75:84)]),
                sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 85:110)],na.rm=TRUE)/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 85:110)],na.rm=TRUE))
            hmd.ax <- nax(cnt.lt[which(cnt.lt$Year==min(y,max(cnt.lt$Year))),],age)
            e0.s <- get.lt(mx=mx.y,age=age,sex="M",hmd.ax)$ex[1]#life expectancy according to stmf data
            e0.h <- get.lt(mx=m.hmd,age=age,sex="M",hmd.ax)$ex[1]#life expectancy according to hmd data
            line <- c(cnt,y,w,e0.s,e0.h)
            mymatrix <- rbind(mymatrix,line)
        }
    }
}

mymatrix <- mymatrix[-1,]
mymatrix.sub <- mymatrix[mymatrix$Year==2019,c(1,3,4)]
colnames(mymatrix.sub)[3] <- "e0.week.2019"
mymatrix2 <- merge(mymatrix,mymatrix.sub,by=c("Country","week"),all.x=TRUE)

mymatrix2$week <- as.numeric(mymatrix2$week)

mymatrix2$Year <- as.factor(mymatrix2$Year)
mymatrix2$e0.week <- as.numeric(mymatrix2$e0.week)
mymatrix2$e0.week.2019 <- as.numeric(mymatrix2$e0.week.2019)

mymatrix2$diff <- mymatrix2$e0.week.2019-mymatrix2$e0.week


library(ggplot2)

# New facet label names for country variable
Country.labs <- c("Austria","Belgium","Bulgaria","Canada","Chile","Croatia","Czechia","Denmark","Estonia","Finland","France","Germany","Greece","Hungary","Iceland","Israel","Italy","Latvia","Lithuania","Luxembourg","Netherlands","New Zealand","Norway","Poland","Portugal","Slovakia","Slovenia","South Korea","Spain","Sweden","Switzerland","England & Wales","Scotland","N. Ireland","USA")
names(Country.labs) <- c("AUT", "BEL", "BGR","CAN","CHL","HRV","CZE","DNK","EST","FIN","FRATNP","DEUTNP","GRC","HUN","ISL","ISR","ITA","LVA","LTU","LUX","NLD","NZL_NP","NOR","POL","PRT","SVK","SVN","KOR","ESP","SWE","CHE","GBRTENW","GBR_SCO","GBR_NIR","USA")
plot1 <- ggplot(data=subset(mymatrix2,Year%in%c("2021")),aes(x=week,y=diff,group=Country,col=Year), show.legend = FALSE)+geom_line(show.legend = FALSE)+ facet_wrap( ~ Country, labeller = labeller(Country=Country.labs)) +ylab("Weekly life expectancy drop")+ geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+scale_color_manual(values=c('Blue'))

##Read Vaccination data from OurWorldInData

vacc <- read.csv("owid-covid-data.csv",sep=",",header=TRUE)[,c(3,4,6,15,37,41,43,46,47,48)]

Countries <- c("Austria","Belgium","Bulgaria","Canada","Chile","Croatia","Czechia","Denmark","Estonia","Finland","France","Germany","Greece","Hungary","Iceland","Israel","Italy","Japan","Latvia","Lithuania","Luxembourg","Netherlands","New Zealand","Norway","Poland","Portugal","Russia","Slovakia","Slovenia","South Korea","Spain","Sweden","Switzerland","Taiwan","United Kingdom","United States")

vacc2 <- vacc[which(vacc$location %in% Countries),]

##Credits to Jonas Schoeley for this function
ISOWeekDate2Date <- function (year, week, weekday = 1, offset = 0) {
  require(ISOweek)
  isoweek_string <-
    paste0(
      year, '-W',
      formatC(
        week+offset,
        flag = '0',
        format = 'd',
        digits = 1
      ),
      '-', weekday
    )
  ISOweek2date(isoweek_string)
}
mymatrix2$date <- ISOWeekDate2Date(mymatrix2$Year,mymatrix2$week)

mymatrix2 <- mymatrix2[mymatrix2$Year%in%c("2020","2021"),]

vacc2$location <- as.factor(vacc2$location)
vacc2$date3 <- as.Date(vacc2$date,format="%Y-%m-%d")+21
vacc2$date <- as.Date(vacc2$date,format="%Y-%m-%d")+21
vacc2$date5 <- as.Date(vacc2$date,format="%Y-%m-%d")+35
colnames(vacc2)[1] <- "Country"
levels(vacc2$Country) <- c( "AUT","BEL","BGR","CAN","CHL","HRV","CZE","DNK","EST","FIN","FRATNP","DEUTNP","GRC", "HUN","ISL","ISR","ITA","JPN","LVA","LTU","LUX","NLD","NZL_NP","NOR","POL","PRT","RUS","SVK","SVN","KOR","ESP","SWE","CHE","TWN","GBRTENW","USA")

library(dplyr)
library(lubridate)

new_matrix <- left_join(vacc2,mymatrix2,by=c("Country","date"))

new_matrix2 <- na.omit(new_matrix)
new_matrix2$month <- month(new_matrix2$date)

library(RColorBrewer)

colrs <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
col <- sample(col_vec, 31)
col <- col_vec[1:31]
    
months.labs <- c("August", "September", "October","November", "December")
names(months.labs) <- c("8", "9", "10","11","12")
plot5 <- ggplot(data=subset(new_matrix2,date>"2021-08-01"&date<="2021-11-30"),aes(x=people_fully_vaccinated_per_hundred,y=diff,col=Country))+ geom_point(show.legend = FALSE)+geom_smooth(data=subset(new_matrix2,date>"2021-08-01"&date<="2021-11-30"&Country!="NZL_NP"&Country!="BGR"),aes(x=people_fully_vaccinated_per_hundred,y=diff),method = "lm",inherit.aes=FALSE)+facet_wrap( ~ month, labeller = labeller(month=months.labs)) +ylab("Weekly life expectancy drop")+xlab("share of fully vaccinated population (4 weeks lagged)")+ scale_color_manual(values=col)

