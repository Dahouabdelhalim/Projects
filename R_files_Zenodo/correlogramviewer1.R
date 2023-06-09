library(ggplot2) # きれいなグラフをプロットするためのライブラリを読み込み
library(signal) # ハイカットフィルタのためのライブラリを読み込み
library(gridExtra) # グラフを縦横に並べるためのライブラリを読み込み
library(ggfortify) # 相互相関関数をプロットして並べるためのライブラリを読み込み
library(reshape2) #

sampling.rate <- 7.5 # サンプリング周波数は7.5 Hz
nyquist.freq <- sampling.rate/2 # ナイキスト周波数は定義からサンプリング周波数の半分
filter.order <- 4 # ハイカットに使うButterworthフィルタは4次のものを使う
ccf.lag.max <- 20*sampling.rate # 長い相互相関の表示は20秒まで
thetheme <- theme_classic(6) # この数字を変えるとフォントサイズが変わる
thepointsize <- 0.1 # この数字を変えると散布図の点のサイズが変わる
maxfignum <- 10 # 一ページあたりの図の最大数
onlyheatmap <- FALSE # ヒートマップ以外の散布図やcross-correlogramを表示しないならTRUE、表示するならFALSE
no.p.value <- FALSE # p値のCSVファイルを出力しないならTRUE
analyze.raw <- FALSE # rawデータもヒートマップを作って検定するならTRUE
np <- 31 # ヒストグラムの階調数
palette <- colorRampPalette(c("blue", "white", "red"))(n=np)
unified.y.scale <- FALSE # 時系列の縦軸の範囲をすべて同一にするならばTRUE、チャネルごとに適宜決めるならFALSE
separated.r.t <- TRUE # 最後のグラフでrestとtaskを分けるならTRUE、分けないならFALSE
offsetcount <- 0 # taskの始めと終わりから除去する長さ


filename <- sprintf("ID0%s.csv", readline("Recording ID? "))
phase = readline("Rest (r), Task (t), Interview (i)? ")
chnum1 = as.numeric(readline("Ch A? (1-40) "))
chnum2 = as.numeric(readline("Ch B? (1-40) "))
high.cut.threshold <- as.numeric(readline("High cut threshold? (originally 1Hz) "))
low.cut.threshold <- as.numeric(readline("Low cut threshold? (originally 0.01Hz) "))
#filename <- "ID016.csv"
#chnum1 <- 1
#chnum2 <- 10
#high.cut.threshold <- 1
#low.cut.threshold <- 0.01

channels <- 40
load.and.split <- function(filename) {
    rawdata <- as.data.frame(read.table(filename, sep=",", dec=".", stringsAsFactors=F))
    offset <- 6
    elinenum <- which(rawdata[[1]]=='Event')
    originalch <- rawdata[elinenum-1,offset:(channels+offset-1)]
    rawdata <- rawdata[(elinenum+1):nrow(rawdata),]
    rawdata <- cbind(rawdata[1], as.data.frame(lapply(rawdata[2], as.numeric)), as.data.frame(lapply(rawdata[offset:(offset+channels-1)], as.numeric)))
    channelname <- gsub("[^a-z0-9]*", "", originalch)
    colnames(rawdata) <- c(c('Event', 'Time'), channelname)
    leftside <- 1:10
    rightside <- c(20,18,19,17,15,16,14,12,13,11)
    rawdata <- rawdata[c(1,2,leftside+2,rightside+2,leftside+20+2,rightside+20+2)]
    return(rawdata)
}

pdf("crosscorrelogram.pdf") # PDFファイルに出力

rawcontents <- load.and.split(filename)

cat("Filtering...\\n")

getnum <- function(x) {
    return(gsub("[^0-9]", "", x))
}

bfilter1 <- butter(filter.order, high.cut.threshold/nyquist.freq, type="low") # Butterworthフィルタを作成
bfilter2 <- butter(filter.order, low.cut.threshold/nyquist.freq, type="high") # Butterworthフィルタを作成
ff <- function(x){return(filtfilt(bfilter2,filtfilt(bfilter1,x)))}
filtered.rawcontents <- data.frame(cbind(rawcontents[1:2], lapply(rawcontents[3:(channels+2)], ff)))

lrstrip <- function(x) {
    return(sub("[[:blank:]]*$", "", sub("^[[:blank:]]*", "", x)))
}

contents <- list()
filtered.contents <- list()
thetime <- list()
name.of.contents <- list()
theraw <- rawcontents
thefiltered <- filtered.rawcontents
nstart <- which(sapply(theraw[1], lrstrip)=="Noise-start")
nend <- which(sapply(theraw[1], lrstrip)=="Noise-end")
if (length(nstart)!=length(nend)) {
    cat("Noise-start and Noise-end do not balance.\\n")
    stopifnot(FALSE)
}
if (length(nstart)>0) {
    nsn <- cbind(nstart,nend)
    skip <- 0
    for (i in 1:nrow(nsn)) {
        theraw <- theraw[c(1:(nsn[i,1]-skip-1), (nsn[i,2]-skip+1):nrow(theraw)),]
        skip <- skip-(nsn[i,2]-nsn[i,1]+1)
    }
}
rstart <- which(sapply(theraw[1], lrstrip)=="Rest-start")
rend <- which(sapply(theraw[1], lrstrip)=="Rest-end")
tstart <- which(sapply(theraw[1], lrstrip)=="Task-start")+offsetcount
tend <- which(sapply(theraw[1], lrstrip)=="Task-end")-offsetcount
istart <- which(sapply(theraw[1], lrstrip)=="Interview-start")+offsetcount
iend <- which(sapply(theraw[1], lrstrip)=="Interview-end")-offsetcount
cut.and.push <- function(start, end, r1, r2, clth, lr) {
    ind <- length(contents)+1
    contents[[ind]] <<- data.frame(theraw[start:end,(r1+2):(r2+2)])
    filtered.contents[[ind]] <<- data.frame(thefiltered[start:end,(r1+2):(r2+2)])
    thetime[[ind]] <<- theraw[start:end,2]
    name.of.contents[ind] <<- paste(clth, lr, sep="")
}
cut.and.push(rstart, rend, 1, 10, "rCl", "L")
cut.and.push(rstart, rend, 11, 20, "rCl", "R")
cut.and.push(tstart, tend, 1, 10, "tCl", "L")
cut.and.push(tstart, tend, 11, 20, "tCl", "R")
cut.and.push(istart, iend, 1, 10, "iCl", "L")
cut.and.push(istart, iend, 11, 20, "iCl", "R")

cut.and.push(rstart, rend, 21, 30, "rTh", "L")
cut.and.push(rstart, rend, 31, 40, "rTh", "R")
cut.and.push(tstart, tend, 21, 30, "tTh", "L")
cut.and.push(tstart, tend, 31, 40, "tTh", "R")
cut.and.push(istart, iend, 21, 30, "iTh", "L")
cut.and.push(istart, iend, 31, 40, "iTh", "R")

names(contents) <- name.of.contents
names(filtered.contents) <- name.of.contents
names(thetime) <- name.of.contents

if (unified.y.scale) {
    ymax <- max(as.data.frame(lapply(contents, max))) # 読み込んだ値の最大値
    ymin <- min(as.data.frame(lapply(contents, min))) # 読み込んだ値の最小値
    yrange = ylim(ymin, ymax)
} else {
    yrange = c()
}

for (name1 in name.of.contents) {
    df <- cbind(thetime[[name1]], filtered.contents[[name1]])
    colnames(df) <- c("Time", colnames(filtered.contents[[name1]]))
    write.csv(df, paste("filtered_timeseries_", name1, ".csv", sep=""), row.names=FALSE)
}

the.name <- function(x) {
    return(substr(x, 2, nchar(x)-1))
}

is.self.simul <- function(x, y) {
    return(substr(x, 1, nchar(x)-1)==substr(y, 1, nchar(y)-1))
}

is.simul <- function(x, y) {
    return((substr(x, 1, 1)==substr(y, 1, 1)) && (gsub("[^0-9]*", "", x)==gsub("[^0-9]*", "", y)))
}

hemisphere <- function(x) {
    return(substr(x, nchar(x), nchar(x)))
}

trstate <- function(x) {
    return(substr(x, 1, 1))
}

fillzero <- function(x) {
    if (nchar(x)==1) {
        return(paste("0", x, sep=""))
    }
    return(x)
}

tile <- function(thelist, n1, n2, maxfignum) {
    N1 = ceiling(n1/maxfignum)
    N2 = ceiling(n2/maxfignum)
    for (i1 in 1:N1) {
        for (i2 in 1:N2) {
            newlist <- list()
            max1 <- 0
            max2 <- 0
            for (j1 in 1:maxfignum) {
                for (j2 in 1:maxfignum) {
                    pos1 <- (i1-1)*maxfignum+j1-1
                    pos2 <- (i2-1)*maxfignum+j2-1
                    if ((pos1<n1) && (pos2<n2)) {
                        newlist <- c(newlist, list(thelist[[pos1*n2+pos2+1]]))
                        max1 <- max(max1, j1)
                        max2 <- max(max2, j2 )
                    }
                }
            }
            do.call(grid.arrange, c(newlist, list(ncol=max2, nrow=max1))) # 散布図を縦max1行、横max2列の中に表示
        }
    }
}

num2name <- function(chnum) {
    if (chnum<=10) {
        return(paste(phase, "ClL", sep=""))
    } else if (chnum<=20) {
        return(paste(phase, "ClR", sep=""))
    } else if (chnum<=30) {
        return(paste(phase, "ThL", sep=""))
    } else {
        return(paste(phase, "ThR", sep=""))
    }
}

name1 <- num2name(chnum1)
name2 <- num2name(chnum2)

dat1 <- filtered.contents[[name1]]
dat2 <- filtered.contents[[name2]]

col1 <- sprintf("ch%d", chnum1)
col2 <- sprintf("ch%d", chnum2)
#dat1[[col1]] <- sin(seq(-40*pi, 40*pi, length=length(dat1[[col1]])))
#dat1[[col2]] <- sin(0.3+seq(-40*pi, 40*pi, length=length(dat1[[col2]])))

timeseries.plots <- c(list(ggplot(data=transform(dat1, sec=thetime[[name1]]), aes_string(x='sec',y=col1))+geom_line()+thetheme+yrange+labs(title=name1)), list(ggplot(data=transform(dat2, sec=thetime[[name2]]), aes_string(x='sec',y=col2))+geom_line()+thetheme+yrange+labs(title=name2))) # リストtimeseries.plotにグラフを追加する
do.call(grid.arrange, c(timeseries.plots, list(ncol=1, nrow=2)))

ccf.plots1 <- list() # 短い相互相関のリストを初期化
ccf.plots2 <- list() # 長い相互相関のリストを初期化
ts1 <- ts(data=dat1[[col1]], frequency=sampling.rate) # col1のチャネルから時系列を作成
ts2 <- ts(data=dat2[[col2]], frequency=sampling.rate) # col2のチャネルから時系列を作成
thetitle <- sprintf("%s\\n%s", col1, col2) # グラフのタイトルを作成
ccf.plots <- c(list(autoplot(ccf(ts1, ts2, plot=FALSE), xlab="Lag (sec)", bcol=NULL)+theme_classic(18)+labs(title=thetitle)), list(autoplot(ccf(ts1, ts2, lag.max=ccf.lag.max, plot=FALSE), xlab="Lag (sec)", bcol=NULL)+theme_classic(18)+labs(title=thetitle))) # ラグの最大値を変えて相互相関の図を作成する
tile(ccf.plots, 2, 1, 2)

datos <- ts.union(ts1, ts2)
                                        #sp <- spectrum(datos, xlim=c(0.01,0.1))
par(cex=1.2, cex.axis=1.2, cex.lab=1.2, cex.main=1.2, cex.sub=1.2)
sp <- spectrum(datos, xlab="Frequency (Hz)", ylab="Power (dB)", main="", log="dB", ci=-1)
abspower <- sqrt((sp$spec**2)%*%matrix(c(1,1),2,1))
cat(sprintf("\\n The peak frequency of cross spectral density is %f Hz.\\nFigures are found in crosscorrelogram.pdf.", sp$freq[which.max(abspower)]))

dev.off() # PDFとしてすべて出力
