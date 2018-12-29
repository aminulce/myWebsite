	const fadeElement = basicScroll.create({
		elem: document.querySelector('.fadeElement'),
		from: 'bottom-middle',
		to: 'top-top',
		props: {
			'--o': {
				from: 0.99,
				to: .01
			}
		}
	})
	fadeElement.start()