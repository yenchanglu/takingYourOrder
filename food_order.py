import pyttsx3

list_food = []
list_drink = []

list_item_price = [0] * 100

class ORDER:
	def reset():
		global list_drink, list_food, list_item_order, list_item_price
		list_item_order = [0] * 100

	def welcome():
		bot = Speaking()
		bot.speak("Hello How are you today")
		bot.speak("What would you like to have")
		ORDER.process_order()

	def process_order():
		while True:
			bot = Speaking()
			print("*" * 26 + " ORDER FOOD & DRINK " + "*" * 26)
			print(" |NO| |FOOD NAME|         |PRICE|   |  |NO| |DRINK NAME|        |PRICE|")
			i = 0
			while i < len(list_food) or i < len(list_drink):
				var_space = 1
				if i <= 8:
					var_space = 2

				if i < len(list_food):
					food = " (" + str(i + 1) + ")" + " " * var_space + str(list_food[i]) + "  | "
				else:
					food = " " * 36 + "| "
				if i < len(list_drink):
					drink = "(" + str(41 + i) + ")" + " " + str(list_drink[i])
				else:
					drink = ""
				print(food, drink)
				i += 1

			print("\n (P) PAYMENT           (M) MAIN MENU           (C) CHANGE ORDER          (E) EXIT\n" + "_" * 72)
			input_ops = input("Please Select Your Operation: ").upper()
			bot.speak("Please Select Your Operation")

			if (input_ops == 'P'):
				print("\n" * 10)
				ORDER.bill_please()
				break
			if (input_ops == 'M'):
				print("\n" * 10)
				ORDER.process_order()
				break
			if (input_ops == 'C'):
				ORDER.modify_order()
				break
			if (input_ops == 'E'):
				ORDER.cancel_order()
				break

			try:
				int(input_ops)
				if ((int(input_ops) <= len(list_food) and int(input_ops) > 0) or (int(input_ops) <= len(list_drink) + 40 and int(input_ops) > 40)):
					try:
						print("\n" + "_" * 72 + "\n" + str(list_food[int(input_ops) - 1]))
					except:
						pass
					try:
						print("\n" + "_" * 72 + "\n" + str(list_drink[int(input_ops) - 41]))
					except:
						pass

					input_qty = input("How Many You Want to Order? ").upper()
					bot.speak("How Many You Want to Order")
					if int(input_qty) != 0:
						if int(input_qty) > 0:
							list_item_order[int(input_ops) - 1] += int(input_qty)
							print("\n" * 10)
							print("Successfully Ordered!")
							ORDER.process_order()
							break
						else:
							if (list_item_order[int(input_ops) - 1]) <= 0:
								print("\n" * 10 + "ERROR: Invalid Input (" + str(input_qty) + "). Try again!")
							else:
								if ( int(input_qty) + list_item_order[int(input_ops) - 1] >= 0):
									list_item_order[int(input_ops) - 1] += int(input_qty)
									print("\n" * 10)
									print("Successfully Ordered!")
									ORDER.process_order()
									break
								else:
									print("\n" * 10 + "ERROR: Invalid Input (" + str(input_qty) + "). Try again!")
					else:
						print("\n" * 10 + "ERROR: Invalid Input (" + str(input_qty) + "). Try again!")
			except:
				print("\n" * 10 + "ERROR: Invalid Input (" + str(input_ops) + "). Try again!")

	def menu_reader():
		file_food = open('menu/list_food.fsd', 'r')
		for i in file_food:
			list_food.append(str(i.strip()))
		file_food.close()

		file_drink = open('menu/list_drink.fsd', 'r')
		for i in file_drink:
			list_drink.append(str(i.strip()))
		file_drink.close()

		i = 0
		while i <= (len(list_food) - 1):
			if 'US' in list_food[i]:
				list_food[i] = str(list_food[i][:list_food[i].index('US') - 1]) + ' ' * (20 - (list_food[i].index('US') - 1)) + str(list_food[i][list_food[i].index('US'):])
			i += 1

		i = 0
		while i <= (len(list_drink) - 1):
			if 'US' in list_drink[i]:
				list_drink[i] = str(list_drink[i][:list_drink[i].index('US') - 1]) + ' ' * (20 - (list_drink[i].index('US') - 1)) + str(list_drink[i][list_drink[i].index('US'):])
			i += 1

	def price_reader():
		global list_food, list_drink
		list_food = sorted(list_food)
		list_drink = sorted(list_drink)

		i = 0
		while i < len(list_food):
			list_item_price[i] = float(list_food[i][int(list_food[i].index("US") + 3):])
			i += 1

		i = 0
		while i < len(list_drink):
			list_item_price[40 + i] = float(list_drink[i][int(list_drink[i].index("US") + 3):])
			i += 1

	def modify_order():
		print("\n" * 10)
		bot = Speaking()
		bot.speak("Here is your order list")
		print("*" * 30 + " ORDER LIST " + "*" * 30 + "\n")
		total_price = 0
		i = 0
		while i < len(list_item_order):
			if(list_item_order[i] != 0):
				if (i >= 0) and (i < 40):
					print(" " * 17 + " (" + str(i + 1) + ")" + str(list_food[i]) + "  x  " + str(list_item_order[i]))
					total_price += list_item_price[i] * list_item_order[i]
				if (i >= 40) and (i < 80):
					print(" " * 17 + " (" + str(i + 1) + ")" + str(list_drink[i - 40]) + "   x  " + str(list_item_order[i]))
					total_price += list_item_price[i] * list_item_order[i]
				i += 1
			else:
				i += 1
		print(" " * 17 + "_" * 35 + "\n" + " " * 17 + "TOTAL PRICES:       US " + str(round(total_price, 2)) + "\n")
		ORDER.process_order()

	def cancel_order():
		bot = Speaking()
		bot.speak("Thank you")
		print("*" * 32 + " THANK YOU " + "*" * 31 + "\n")

	def bill_please():
		while True:
			bot = Speaking()
			bot.speak("Here is your bill")
			print("*" * 32 + " PAYMENT " + "*" * 33 + "\n")
			total_price = 0
			i = 0
			while i < len(list_item_order):
				if(list_item_order[i] != 0):
					if (i >= 0) and (i < 40):
						print(" " * 17 + str(list_food[i]) + "  x  " + str(list_item_order[i]))
						total_price += list_item_price[i] * list_item_order[i]
					if (i >= 40) and (i < 80):
						print(" " * 17 + str(list_drink[i - 40]) + "   x  " + str(list_item_order[i]))
						total_price += list_item_price[i] * list_item_order[i]
					i += 1
				else:
					i += 1
			print(" " * 17 + "_" * 35 + "\n" + " " * 17 + "TOTAL PRICES:       US " + str(round(total_price, 2)))
			print("\n (P) PAY           (M) MAIN MENU           (C) CHANGE ORDER          (E) EXIT\n" + "_" * 72)
			bot.speak("Please select your operation")
			input_ops = str(input("Please Select Your Operation: ")).upper()
			if (input_ops == 'P'):
				print("\n" * 10)
				if (total_price > 0):
					bot.speak("Successfully paid")
					print("Successfully Paid!")
				else:
					bot.speak("Thank you")
					print("*" * 32 + " THANK YOU " + "*" * 31 + "\n")

				break

			elif (input_ops == 'M'):
				print("\n" * 10)
				ORDER.process_order()
				break
			elif (input_ops == 'C'):
				print("\n" * 10)
				ORDER.modify_order()
				break
			elif ('E' in input_ops) or ('e' in input_ops):
				ORDER.cancel_order()
				break
			else:
				print("\n" * 10 + "ERROR: Invalid Input (" + str(input_ops) + "). Try again!")

class Speaking:
	def speak(self, audio):
		engine = pyttsx3.init()
		engine.say(audio)
		engine.runAndWait()

if __name__ == '__main__':
	ORDER()
	ORDER.reset()
	ORDER.menu_reader()
	ORDER.price_reader()
	ORDER.welcome()