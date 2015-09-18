
# source("COZIR.R", print.eval = TRUE)

system("MODE COM8:9600, N, 8, 1")
cozir = file("COM8", open = 'r+')
isOpen(cozir)

# Set operting mode to STOP
write("K 0\r\n", cozir)
Sys.sleep(0.1)
readLines(cozir, n = -1)

# Return firmware version and sensor serial number

write("Y\r\n", cozir)
Sys.sleep(0.1)
s <- readLines(cozir, n = -1)
print(s)
for(i in 1:20){
  write("T\r\n", cozir) # Return most recent temperature measure
  Sys.sleep(0.1)
  t <- scan(cozir, n = 3, what = "character")
  temperature <- as.numeric(t[2])
  print(temperature)
  write("H\r\n", cozir) # Return most recent humidity measure
  Sys.sleep(0.1)
  h <- scan(cozir, n = 3, what = "character")
  humidity <- as.numeric(h[2])
  print(humidity)
  write("Z\r\n", cozir) # Return most recent CO2 measure
  Sys.sleep(0.1)
  c <- scan(cozir, n = 3, what = "character")
  co2 <- as.numeric(c[2])
  print(co2)
  datarow <- c(temperature, humidity, co2)
  if (i == 1){
    cozirData <- datarow
  } else {
    cozirData <- rbind(cozirData, datarow)
  }
}

close(cozir)